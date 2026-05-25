#!/usr/bin/env python3
"""
bisulfite_amplicon_analysis.py
Demirer Lab — Caltech CCE

CyMATE-style bisulfite amplicon methylation analysis for Oxford Nanopore reads.

Performs bisulfite-aware semi-global alignment of nanopore amplicon reads to a
genomic plus-strand reference, calls per-site CpG / CHG / CHH methylation, and
optionally computes epiallele patterns. Outputs:
  • CyMATE-compatible TSV (position, context, C/T read counts, %mC)
  • Publication-quality per-site bar chart (PNG + SVG; auto-split if amplicon is long)
  • Composite figure with epiallele heatmap above bar chart (if EPIALLELE = True)

Dependencies: numpy, matplotlib (standard scientific Python stack)

Usage:
    Edit the CONFIG block below, then:
        python3 bisulfite_amplicon_analysis.py
"""

import csv, os, random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG — edit this block for each amplicon / sample
# ══════════════════════════════════════════════════════════════════════════════

# Genomic plus-strand reference sequence spanning both primers (5′→3′).
# Must start with the forward primer sequence and end with the reverse primer
# sequence (on the plus strand, i.e. the RC of the actual reverse primer oligo).
# Spaces are stripped automatically — use them freely for readability.
TEMPLATE = (
    'AAGGCCGAAGAGGGGAAAGGTTCCATGTGAACGGCACTTGCACATGGGTTAGTCGATCCT'
    'AAGAGTCGGGGGAAACCCGTCTGATAGCGCTTAAGCGCGAACTTCGAAAGGGGATCCGGT'
    'TAAAATTCCGGAACCGGGACGTGGCGGTTGACGGCAACGTTAGGGAGTCCGGAGACGTCGG'
    'CGGGGGCCTCGGGAAGAGTTATCTTTTCTGTTTAACAGCCTGCCCACCCTGGAAACGGCTC'
    'AGCCGGAGGTAGGGTCCAGCGGCTGGAAGAGCACCGCACGTCGCGTGGTGTCCGGTGCGCC'
    'CCCGGCGGCCCTTGAAAATCCGGAGGACCGAGTGCCGCTCACGCCCGGTCGTACTCATAAC'
    'CGCATCAGGTCTCCAAGGTGAACAGCCTCTGGTCGATGGAACAATGTAGGCAAGGGAAGT'
).replace(' ', '')

# Primer sequences (IUPAC degenerate bases supported; Y=C/T, R=A/G, etc.).
# FWD_PRIMER: 5′→3′ sequence as synthesised (binds plus strand).
# REV_PRIMER: 5′→3′ sequence as synthesised (binds minus strand).
FWD_PRIMER = 'AAGGYYGAAGAGGGGAAAGG'      # RE0708-rDNAmeth4F  (20 bp)
REV_PRIMER = 'ACTTCCCTTRCCTACATTRTTCCA'  # RE0709             (24 bp)

# Input FASTQ and output settings
FASTQ  = 'path/to/sample.fastq'
OUTDIR = 'output'
LABEL  = 'sample_amplicon'   # used as filename prefix and figure title

# Read-length filter (bp of primer-trimmed sequence).
# Set to ~85–90 % of expected amplicon length to exclude truncated reads.
# Increase for longer amplicons with bimodal length distributions.
MIN_READ_LEN = 50   # set higher (e.g. 300, 700, 900) to enforce full-length reads

# Epiallele analysis — set True only for shorter amplicons (≲ 600 bp) where
# individual reads span the full amplicon and epiallele patterns are meaningful.
EPIALLELE        = True
EPIALLELE_MINLEN = 350   # minimum trimmed read length for epiallele inclusion
EPIALLELE_TOP_N  = 10    # number of top epiallele patterns to display
MAX_NOCALL_FRAC  = 0.20  # exclude reads with > this fraction of no-calls at C sites

# Figure y-axis maximum (%). Set None for auto.
YMAX = None

# ══════════════════════════════════════════════════════════════════════════════
# IUPAC / reverse-complement utilities
# ══════════════════════════════════════════════════════════════════════════════

# Bisulfite-aware IUPAC: C in primer matches C or T in bisulfite-converted read
IUPAC_BS = {
    'A': 'A',   'C': 'CT',  'G': 'G',   'T': 'T',
    'R': 'AG',  'Y': 'CT',  'M': 'AC',  'K': 'GT',
    'S': 'CG',  'W': 'AT',  'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}

_DNA_RC   = str.maketrans('ACGT', 'TGCA')
_DEGEN_RC = str.maketrans('ACGTRYMKSWBDHVacgtrymkswbdhv',
                           'TGCAYRKMSWVHDBtgcayrkmswvhdb')

def rc(seq):
    """Reverse-complement of an unambiguous DNA sequence."""
    return seq.translate(_DNA_RC)[::-1]

def rc_degen(seq):
    """Reverse-complement of a degenerate (IUPAC) primer sequence."""
    return seq.translate(_DEGEN_RC)[::-1]

# ══════════════════════════════════════════════════════════════════════════════
# Amplicon geometry (derived from TEMPLATE + primers)
# ══════════════════════════════════════════════════════════════════════════════

TEMPLATE       = TEMPLATE.replace(' ', '')
REV_PRIMER_RC  = rc_degen(REV_PRIMER)
AMP_START      = len(FWD_PRIMER)
AMP_END        = len(TEMPLATE) - len(REV_PRIMER)
AMPLICON       = TEMPLATE[AMP_START:AMP_END]   # inter-primer sequence

# ══════════════════════════════════════════════════════════════════════════════
# I/O
# ══════════════════════════════════════════════════════════════════════════════

def parse_fastq(path):
    """Return list of upper-case read sequences from a FASTQ file."""
    reads = []
    with open(path) as fh:
        lines = fh.readlines()
    for i in range(0, len(lines) - 3, 4):
        seq = lines[i + 1].strip().upper()
        if seq:
            reads.append(seq)
    return reads

# ══════════════════════════════════════════════════════════════════════════════
# Strand detection (bisulfite C/G ratio heuristic)
# ══════════════════════════════════════════════════════════════════════════════

def bisulfite_strand(seq):
    """
    Infer bisulfite strand from C and G frequency.
    Plus strand (+): most cytosines converted → C-depleted, G-rich.
    Minus strand (−): RC'd plus strand → G-depleted, C-rich.
    Returns '+', '-', or '?' if ambiguous.
    """
    n = len(seq)
    if n == 0:
        return '?'
    c = seq.count('C') / n
    g = seq.count('G') / n
    if c < 0.15 and c + 0.08 < g:
        return '+'
    if g < 0.15 and g + 0.08 < c:
        return '-'
    return '?'

# ══════════════════════════════════════════════════════════════════════════════
# Primer detection (degenerate, bisulfite-aware, ≤ 3 mismatches)
# ══════════════════════════════════════════════════════════════════════════════

def _deg_mm(primer, segment):
    """Count mismatches between a degenerate primer and a read segment."""
    mm = 0
    for p, q in zip(primer, segment):
        allowed = IUPAC_BS.get(p, p)
        if 'C' in allowed and q == 'T':   # bisulfite C→T is a valid match
            continue
        if q not in allowed:
            mm += 1
    return mm

def find_fwd(read, max_offset=15, max_mm=3):
    """Return start index of FWD primer in read, or None."""
    plen = len(FWD_PRIMER)
    for off in range(min(max_offset + 1, max(0, len(read) - plen + 1))):
        if _deg_mm(FWD_PRIMER, read[off:off + plen]) <= max_mm:
            return off
    return None

def find_rev_rc(read, max_offset=15, max_mm=3):
    """Return start index of REV primer RC in read (at 3′ end), or None."""
    plen = len(REV_PRIMER_RC)
    for off in range(min(max_offset + 1, max(0, len(read) - plen + 1))):
        pos = len(read) - plen - off
        if pos >= 0 and _deg_mm(REV_PRIMER_RC, read[pos:pos + plen]) <= max_mm:
            return pos
    return None

# ══════════════════════════════════════════════════════════════════════════════
# Bisulfite-aware semi-global alignment (read inside reference)
# ══════════════════════════════════════════════════════════════════════════════

def align_bisulfite(ref, read, MATCH=2, MISMATCH=-2, GAP=-2):
    """
    Semi-global dynamic-programming alignment: the read is fully consumed
    (global in the read dimension) but the reference may have free leading
    and trailing gaps (local in the reference dimension).  This allows short
    or partial reads to align to the correct sub-region of the amplicon.

    Bisulfite scoring: C (ref) vs T (read) is treated as a match, reflecting
    unmethylated cytosine conversion.

    Returns
    -------
    aref      : str  aligned reference string (gaps as '-')
    aread     : str  aligned read string (gaps as '-')
    ref_start : int  0-based start position in ref
    ref_end   : int  0-based exclusive end position in ref
    """
    m, n = len(ref), len(read)
    r_b = np.frombuffer(ref.encode('ascii'),  dtype=np.uint8)
    q_b = np.frombuffer(read.encode('ascii'), dtype=np.uint8)

    exact = (r_b[:, None] == q_b[None, :])
    bs    = (r_b[:, None] == ord('C')) & (q_b[None, :] == ord('T'))
    S     = np.where(exact | bs, float(MATCH), float(MISMATCH)).astype(np.float32)

    dp = np.empty((m + 1, n + 1), dtype=np.float32)
    dp[:, 0] = 0.0                         # free leading reference gaps
    dp[0, :] = np.arange(n + 1) * GAP     # read must be consumed from start

    j_idx = np.arange(1, n + 1, dtype=np.float32)
    k_idx = np.arange(n,        dtype=np.float32)

    for i in range(1, m + 1):
        diag    = dp[i - 1, :-1] + S[i - 1]
        up      = dp[i - 1, 1:]  + GAP
        best_du = np.maximum(diag, up)
        v       = best_du + k_idx * (-GAP)
        rmax    = np.maximum.accumulate(v)
        dp[i, 1:] = np.maximum(rmax + (j_idx - 1) * GAP,
                                dp[i, 0] + j_idx * GAP)

    best_i = int(np.argmax(dp[:, n]))   # best reference end position

    # Traceback
    aref, aread = [], []
    i, j = best_i, n
    EPS  = 0.5
    while i > 0 or j > 0:
        if i == 0:
            aref.append('-'); aread.append(read[j - 1]); j -= 1
            continue
        if j == 0:
            break
        r, q = ref[i - 1], read[j - 1]
        sc   = MATCH if (r == q or (r == 'C' and q == 'T')) else MISMATCH
        if abs(dp[i, j] - (dp[i - 1, j - 1] + sc)) < EPS:
            aref.append(r);          aread.append(q);       i -= 1; j -= 1
        elif abs(dp[i, j] - (dp[i - 1, j] + GAP)) < EPS:
            aref.append(ref[i - 1]); aread.append('-');     i -= 1
        else:
            aref.append('-');        aread.append(read[j - 1]); j -= 1

    aref  = list(reversed(aref))
    aread = list(reversed(aread))
    ref_end   = best_i
    ref_start = ref_end - sum(1 for c in aref if c != '-')
    return ''.join(aref), ''.join(aread), ref_start, ref_end

# ══════════════════════════════════════════════════════════════════════════════
# Methylation context
# ══════════════════════════════════════════════════════════════════════════════

def get_context(seq, pos):
    """Return 'CpG', 'CHG', or 'CHH' for a cytosine at pos in seq."""
    if pos + 1 < len(seq) and seq[pos + 1] == 'G':
        return 'CpG'
    if pos + 2 < len(seq) and seq[pos + 2] == 'G':
        return 'CHG'
    return 'CHH'

# ══════════════════════════════════════════════════════════════════════════════
# Main analysis
# ══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUTDIR, exist_ok=True)

print(f'Template : {len(TEMPLATE)} bp')
print(f'Amplicon : {len(AMPLICON)} bp  (positions {AMP_START}–{AMP_END})')

reads_raw = parse_fastq(FASTQ)
print(f'\nTotal reads  : {len(reads_raw)}')

# Strand orientation
oriented, strand_ctr = [], Counter()
for seq in reads_raw:
    st = bisulfite_strand(seq)
    strand_ctr[st] += 1
    if   st == '+': oriented.append(seq)
    elif st == '-': oriented.append(rc(seq))

print(f'Strand calls : {dict(strand_ctr)}')
print(f'Oriented     : {len(oriented)}')

amp_len   = len(AMPLICON)
c_pos     = [i for i, b in enumerate(AMPLICON) if b == 'C']
c_pos_set = set(c_pos)
contexts  = {i: get_context(AMPLICON, i) for i in c_pos}

C_counts   = np.zeros(amp_len, dtype=np.int32)
T_counts   = np.zeros(amp_len, dtype=np.int32)
read_meths = []   # list of (meth_dict, trimmed_len) — used for epiallele analysis

n_primer = n_aligned = n_skipped = 0

for seq in oriented:
    fwd_pos = find_fwd(seq)
    if fwd_pos is None:
        n_skipped += 1
        continue
    n_primer += 1

    trim_start = fwd_pos + len(FWD_PRIMER)
    rev_pos    = find_rev_rc(seq)
    trim_end   = rev_pos if rev_pos is not None else len(seq)
    trimmed    = seq[trim_start:trim_end]

    if len(trimmed) < MIN_READ_LEN:
        n_skipped += 1
        continue

    aref, aread, ref_start, _ = align_bisulfite(AMPLICON, trimmed)

    meth    = {}
    ref_pos = ref_start
    for r, q in zip(aref, aread):
        if r != '-':
            if ref_pos in c_pos_set:
                if   q == 'C': C_counts[ref_pos] += 1; meth[ref_pos] = 'M'
                elif q == 'T': T_counts[ref_pos] += 1; meth[ref_pos] = 'U'
                else:          meth[ref_pos] = '?'
            ref_pos += 1

    read_meths.append((meth, len(trimmed)))
    n_aligned += 1

print(f'Primer found : {n_primer}')
print(f'Aligned      : {n_aligned}')
print(f'Skipped      : {n_skipped}')

# ── Per-site results ──────────────────────────────────────────────────────────

results = []
for ap in c_pos:
    tp    = ap + AMP_START        # 0-based template position
    n_c   = int(C_counts[ap])
    n_t   = int(T_counts[ap])
    total = n_c + n_t
    meth  = 100.0 * n_c / total if total > 0 else None
    results.append({
        'pos':            tp,
        'context':        contexts[ap],
        'C_reads':        n_c,
        'T_reads':        n_t,
        'total_reads':    total,
        'methylation_pct': meth,
        'covered':        total > 0,
    })

print('\nPer-context methylation summary:')
for ctx in ['CpG', 'CHG', 'CHH']:
    rows = [r for r in results if r['context'] == ctx and r['covered']]
    if rows:
        vals = [r['methylation_pct'] for r in rows]
        print(f'  {ctx:3s} : {len(rows):3d} sites  '
              f'mean={np.mean(vals):.1f}%  '
              f'range={min(vals):.1f}–{max(vals):.1f}%')

# ── TSV output ────────────────────────────────────────────────────────────────

tsv_path = os.path.join(OUTDIR, f'{LABEL}_methylation_cymate.tsv')
with open(tsv_path, 'w', newline='') as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=['pos', 'context', 'C_reads', 'T_reads',
                    'total_reads', 'methylation_pct', 'covered'],
        delimiter='\t',
    )
    w.writeheader()
    for r in results:
        row = dict(r)
        pct = row['methylation_pct']
        row['methylation_pct'] = '' if pct is None else f'{pct:.4f}'
        w.writerow(row)
print(f'\nSaved TSV : {tsv_path}')

# ══════════════════════════════════════════════════════════════════════════════
# Epiallele analysis (optional)
# ══════════════════════════════════════════════════════════════════════════════

epi_patterns, epi_ctr, n_epi = [], Counter(), 0

if EPIALLELE:
    n_c_pos = len(c_pos)
    for meth, rlen in read_meths:
        if rlen < EPIALLELE_MINLEN:
            continue
        pat = tuple(meth.get(p, '?') for p in c_pos)
        if pat.count('?') <= n_c_pos * MAX_NOCALL_FRAC:
            epi_patterns.append(pat)

    epi_ctr  = Counter(epi_patterns)
    n_epi    = len(epi_patterns)
    n_unique = len(epi_ctr)
    print(f'\nEpiallele reads    : {n_epi}')
    print(f'Unique epialleles  : {n_unique}')

# ══════════════════════════════════════════════════════════════════════════════
# Figure helpers
# ══════════════════════════════════════════════════════════════════════════════

BAR_CLR = {'CpG': '#C1440E', 'CHG': '#3AAFA9', 'CHH': '#9B72AA'}

legend_patches = [
    mpatches.Patch(color=BAR_CLR['CpG'], label='CpG'),
    mpatches.Patch(color=BAR_CLR['CHG'], label='CHG'),
    mpatches.Patch(color=BAR_CLR['CHH'], label='CHH'),
]

all_pos    = sorted(r['pos'] for r in results if r['covered'])
pos_to_idx = {p: i for i, p in enumerate(all_pos)}
ctx_data   = {ctx: {} for ctx in ['CpG', 'CHG', 'CHH']}
for r in results:
    if r['covered']:
        ctx_data[r['context']][r['pos']] = r['methylation_pct']

# Auto y-axis ceiling
if YMAX is None:
    all_vals = [v for d in ctx_data.values() for v in d.values()]
    YMAX = max(int(np.ceil(max(all_vals) / 5) * 5 * 1.15), 20) if all_vals else 60


def _draw_bar_panel(ax, positions, data_dict, xlabel=False, legend=False, title=''):
    """Draw a single per-site mC bar panel onto ax."""
    for i, p in enumerate(positions):
        ctx = next((c for c in ['CpG', 'CHG', 'CHH'] if p in data_dict[c]), None)
        if ctx is None:
            continue
        ax.bar(i, data_dict[ctx][p], width=0.80, color=BAR_CLR[ctx], zorder=2)

    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([str(p - AMP_START) for p in positions],
                        rotation=45, ha='right', fontsize=11)
    for lbl in ax.get_xticklabels():
        lbl.set_fontweight('bold'); lbl.set_fontfamily('Arial')
    ax.tick_params(axis='x', length=6, width=1.2, direction='out', pad=4)
    ax.tick_params(axis='y', labelsize=20)
    for lbl in ax.get_yticklabels():
        lbl.set_fontweight('bold'); lbl.set_fontfamily('Arial')
    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylim(0, YMAX)
    yticks = [0] + [t for t in [15, 30, 45, 60, 75, 100] if t <= YMAX]
    ax.set_yticks(yticks)
    ax.set_ylabel('% mC', fontsize=22, fontweight='bold', labelpad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if title:
        ax.text(0.005, 0.95, title, transform=ax.transAxes,
                ha='left', va='top', fontsize=20, fontweight='bold',
                color='#333333', fontfamily='Arial')
    if xlabel:
        ax.set_xlabel('Position (amplicon bp)', fontsize=22,
                      fontweight='bold', labelpad=16)
    if legend:
        ax.legend(handles=legend_patches, fontsize=17,
                  loc='upper right', framealpha=0.9)

plt.rcParams['font.family'] = 'Gill Sans'

# ══════════════════════════════════════════════════════════════════════════════
# Figure A — per-site bar chart
# Automatically split into two panels when there are > 150 covered positions.
# ══════════════════════════════════════════════════════════════════════════════

SPLIT_THRESHOLD = 150
n_panels = 2 if len(all_pos) > SPLIT_THRESHOLD else 1

if n_panels == 1:
    fig, ax = plt.subplots(1, 1, figsize=(32, 5), dpi=200)
    _draw_bar_panel(ax, all_pos, ctx_data,
                    xlabel=True, legend=True,
                    title=f'Per-site cytosine methylation — {LABEL}')
else:
    mid      = len(all_pos) // 2
    halves   = [all_pos[:mid], all_pos[mid:]]
    fig, axes = plt.subplots(2, 1, figsize=(32, 9), dpi=200)
    _draw_bar_panel(axes[0], halves[0], ctx_data,
                    legend=True,
                    title=f'Per-site cytosine methylation — {LABEL}  (part 1)')
    _draw_bar_panel(axes[1], halves[1], ctx_data,
                    xlabel=True,
                    title=f'Per-site cytosine methylation — {LABEL}  (part 2)')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.60)

for ext in ('png', 'svg'):
    p = os.path.join(OUTDIR, f'{LABEL}_mC_barchart.{ext}')
    plt.savefig(p, dpi=400, bbox_inches='tight')
    print(f'Saved : {p}')
plt.close()

# ══════════════════════════════════════════════════════════════════════════════
# Figure B — composite (epiallele heatmap + bar chart)
# Only produced when EPIALLELE = True and enough patterns were found.
# ══════════════════════════════════════════════════════════════════════════════

if EPIALLELE and n_epi >= 2:
    from matplotlib.transforms import blended_transform_factory

    n_c_pos     = len(c_pos)
    top_N       = min(EPIALLELE_TOP_N, len(epi_ctr))
    top_patterns = epi_ctr.most_common(top_N)
    freqs        = [cnt / n_epi * 100     for _, cnt in top_patterns]
    counts_top   = [cnt                   for _, cnt in top_patterns]

    # Context color strip above heatmap
    CTX_CLR = {
        'CpG': [0.757, 0.267, 0.055],   # #C1440E
        'CHG': [0.227, 0.686, 0.663],   # #3AAFA9
        'CHH': [0.608, 0.447, 0.667],   # #9B72AA
    }
    ctx_strip = np.array([[CTX_CLR[contexts[c_pos[j]]]
                            for j in range(n_c_pos)]], dtype=np.float32)

    # Binary matrix: M=black, U=white, ?=light-grey
    binary = np.ones((top_N, n_c_pos, 3), dtype=np.float32)
    for ri, (pat, _) in enumerate(top_patterns):
        for ci, call in enumerate(pat):
            if   call == 'M': binary[ri, ci] = [0.0, 0.0, 0.0]
            elif call == 'U': binary[ri, ci] = [1.0, 1.0, 1.0]
            else:             binary[ri, ci] = [0.88, 0.88, 0.88]

    comp_legend = [
        mpatches.Patch(facecolor='black',                        label='Methylated'),
        mpatches.Patch(facecolor='white', edgecolor='#999999',
                       linewidth=0.8,                            label='Unmethylated'),
        mpatches.Patch(facecolor=[0.88, 0.88, 0.88],             label='No call'),
        mpatches.Patch(color=CTX_CLR['CpG'],                     label='CpG'),
        mpatches.Patch(color=CTX_CLR['CHG'],                     label='CHG'),
        mpatches.Patch(color=CTX_CLR['CHH'],                     label='CHH'),
    ]

    ROW_H    = 0.28
    STRIP_H  = 0.14
    LEG_H    = 0.70
    BAR_H    = 4.00
    SPACER_H = 0.15
    TOTAL_H  = LEG_H + STRIP_H + top_N * ROW_H + SPACER_H + BAR_H + 2.0

    fig = plt.figure(figsize=(32, TOTAL_H), dpi=200)
    gs  = fig.add_gridspec(
        5, 2,
        height_ratios=[LEG_H, STRIP_H, top_N * ROW_H, SPACER_H, BAR_H],
        width_ratios=[n_c_pos, 5],
        hspace=0.0, wspace=0.02,
        left=0.06, right=0.97,
        top=0.97, bottom=0.12,
    )

    ax_leg    = fig.add_subplot(gs[0, 0])
    ax_strip  = fig.add_subplot(gs[1, 0])
    ax_epi    = fig.add_subplot(gs[2, 0])
    ax_spacer = fig.add_subplot(gs[3, 0])
    ax_bar    = fig.add_subplot(gs[4, 0])
    ax_freq   = fig.add_subplot(gs[2, 1])
    ax_spacer.axis('off')
    for r in [0, 1, 3, 4]:
        fig.add_subplot(gs[r, 1]).axis('off')

    # Legend
    ax_leg.axis('off')
    ax_leg.legend(handles=comp_legend, fontsize=21, ncol=6,
                  loc='center left', framealpha=0.9,
                  bbox_to_anchor=(0.0, 0.5))
    ax_leg.text(0.72, 0.5,
                f'Top {top_N} of {len(epi_ctr)} unique epialleles shown',
                transform=ax_leg.transAxes, ha='left', va='center',
                fontsize=24, color='#333333')

    # Context strip
    ax_strip.imshow(ctx_strip, aspect='auto', interpolation='nearest')
    ax_strip.set_xlim(-0.5, n_c_pos - 0.5)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    for sp in ax_strip.spines.values(): sp.set_visible(False)
    ax_strip.set_title(LABEL, fontsize=30, fontweight='bold', pad=12)

    # Epiallele heatmap
    ax_epi.imshow(binary, aspect='auto', interpolation='nearest', origin='lower')
    for y in np.arange(-0.5, top_N, 1):
        ax_epi.axhline(y, color='#666666', linewidth=0.6, zorder=3)
    for x in np.arange(-0.5, n_c_pos, 1):
        ax_epi.axvline(x, color='#888888', linewidth=0.8, zorder=3)
    ax_epi.set_xticks(range(n_c_pos))
    ax_epi.set_xticklabels(['' for _ in range(n_c_pos)])
    ax_epi.tick_params(axis='x', length=0, pad=0)
    ax_epi.set_yticks(range(top_N))
    ax_epi.set_yticklabels([f'#{i + 1}' for i in range(top_N)], fontsize=24)
    ax_epi.set_ylabel('Epiallele rank', fontsize=27)
    ax_epi.set_xlim(-0.5, n_c_pos - 0.5)
    for sp in ax_epi.spines.values(): sp.set_visible(False)

    # Read counts / frequency labels
    ax_freq.axis('off')
    trans = blended_transform_factory(ax_epi.transAxes, ax_epi.transData)
    for yi, (f, c) in enumerate(zip(freqs, counts_top)):
        ax_epi.text(1.003, yi, f'{c}  ({f:.1f}%)',
                    transform=trans, ha='left', va='center',
                    fontsize=22, color='#222222')

    # Per-site bar chart (bottom panel)
    _draw_bar_panel(ax_bar, all_pos, ctx_data,
                    xlabel=True,
                    title='Per-site cytosine methylation')

    for ext in ('png', 'svg'):
        p = os.path.join(OUTDIR, f'{LABEL}_composite.{ext}')
        plt.savefig(p, dpi=400, bbox_inches='tight')
        print(f'Saved : {p}')
    plt.close()

print('\nDone.')
