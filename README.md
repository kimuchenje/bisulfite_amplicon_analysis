# bisulfite_amplicon_analysis

Processes Oxford Nanopore bisulfite amplicon sequencing data (FASTQ) and produces per-site cytosine methylation calls and figures in the style of CyMATE.

Demirer Lab, Caltech Division of Chemistry and Chemical Engineering (CCE)

---

## What it does

1. Parses FASTQ reads and assigns strand orientation based on C:G frequency
2. Locates and trims primers using bisulfite-aware IUPAC mismatch search
3. Aligns each read to the amplicon reference using a bisulfite-aware semi-global algorithm
4. Calls methylation per site in CpG, CHG, and CHH contexts
5. Optionally computes epiallele patterns across reads
6. Outputs a CyMATE-compatible TSV and PNG/SVG figures

## Outputs

- `{LABEL}_methylation_cymate.tsv` — per-site methylation table
- `{LABEL}_mC_barchart.png/svg` — per-site methylation bar chart colored by context
- `{LABEL}_composite.png/svg` — epiallele heatmap + bar chart (when `EPIALLELE = True`)

## Usage

Edit the `Config` block at the top of the script with your FASTQ path, primer sequences, and reference template, then run:

```bash
python bisulfite_amplicon_analysis.py
```

## Dependencies

```bash
pip install numpy matplotlib
```

Python ≥ 3.8

## License

MIT — see [LICENSE](LICENSE)
