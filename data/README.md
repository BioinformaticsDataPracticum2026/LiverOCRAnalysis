# Datasets and Input Resources

This folder contains all datasets used in the pipeline, including ATAC-seq peak files, quality control outputs, and annotation files required for downstream analysis. These datasets form the foundation for studying cross-species conservation of open chromatin regions (OCRs) between human and mouse.

---

## Folder Structure

data/
├── atac_seq_peaks/
├── quality_control/
├── tss_annotations/

---

## ATAC-seq Peaks (`atac_seq_peaks/`)

- `humanIDRConservedPeaks.gz`
- `mouseIDRConservedPeaks.gz`

These files contain high-confidence OCR peaks after IDR filtering. They represent regions of open chromatin and are used as the starting point for the pipeline. These peaks are later mapped across species using HALPER to construct orthologous regulatory regions.

---

## Quality Control (`quality_control/`)

- `human_liver_qc.html`
- `mouse_liver_qc.html`
- `human_ovary_qc.html`
- `mouse_ovary_qc.html`
- `qc_summary.txt`

These files contain quality control reports and summary statistics for all datasets. QC evaluation was used to decide which datasets to use for downstream analysis. Based on this analysis, the liver datasets were selected because they showed stronger reproducibility and overall signal quality, while ovary datasets showed signs of over-fragmentation and lower consistency.

### Data Quality Metrics

| Metric | Ovary (Human rep1) | Ovary (Human rep2) | Ovary (Mouse rep1) | Ovary (Mouse rep2) | Liver (Human rep1) | Liver (Human rep2) | Liver (Mouse rep1) | Liver (Mouse rep2) |
|--------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|
| % Mapped Reads | 97.9 | 98.2 | 98.5 | 98.6 | 98.9 | 98.6 | 97.8 | 97.9 |
| % Properly Paired Reads | 96.2 | 96.6 | 94.4 | 94.6 | 97.9 | 97.5 | 95.7 | 95.6 |
| % Mitochondrial Reads | 4.77 | 4.44 | 4.56 | 3.29 | 34.81 | 14.22 | 0.60 | 0.55 |
| Filtered Read Counts (million) | 48.5 | 138.3 | 13.7 | 36.2 | 47.8 | 101.5 | 54.7 | 54.8 |
| NRF (%) | 87.6 | 94.0 | 89.0 | 85.2 | 88.5 | 90.0 | 93.6 | 94.4 |
| TSS Enrichment Peak | 6.64 | 13.8 | 17.6 | 13.5 | 23.3 | 21.0 | 7.69 | 7.35 |

**Summary:**
- All samples show high mapping rates (>97%)
- Strong pairing rates indicate good sequencing quality
- Human liver shows higher mitochondrial reads (possible bias)
- TSS enrichment is strongest in liver datasets
- NRF values indicate good library complexity overall

---

## TSS Annotations (`tss_annotations/`)

- `gencode.v41.annotation.bed` (Human - hg38)
- `gencode.v41.annotation.gtf.gz`
- `gencode.vM30.annotation.bed` (Mouse - mm10)
- `gencode.vM30.annotation.gtf.gz`

These annotation files are used to define transcription start sites (TSS). They are required for classifying OCRs into promoters and enhancers:

- Promoters: OCRs within ±2 kb of a TSS
- Enhancers: OCRs outside this region

---

## Notes

- All files in this folder are input datasets for the pipeline
- Downstream outputs are written to the `results/` directory
- Dataset selection (liver vs ovary) is based on QC evaluation