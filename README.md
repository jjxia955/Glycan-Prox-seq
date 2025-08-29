# Glycan-Prox-seq
**Integrated profiling of proteins, glycans, protein-specific glycoforms, and mRNA in single cells with glycan proximity sequencing**
<div align="center">
  <img src="GPS_shceme.png" alt="GPS Scheme" width="600"/>
  <br>
  <em>Glycan-Prox-seq (GPS) workflow</em>
</div>

---

This repository contains code and data to **reproduce the main figures and analyses** in the Glycan-Prox-seq study.

Author: Junjie Xia
</div>
Email: jjxia@uchicago.edu

---

## Repository Structure

- `code_reproduce/` — Jupyter notebooks (`.ipynb`) and R scripts (`.R`) for reproducing figures  
  - `Fig1.ipynb` — Bulk GPS assay on WT and GnTI knock-out Expi293F cell lines
  - `Fig2.ipynb` — Single-cell plate-based GPS assay on Jurkat and Raji mixtures  
  - `Fig3_*.ipynb` — 10x GPS assay on human PBMCs, including clustering and annotation, ADT/LDT profiles, NK subset analysis.
  - `Fig4&5_pseudotime_analysis.ipynb` — pseudotime trajectory analysis on T cell populations
  - `Fig6_*.ipynb` — Glycosylation remodeling under external perturbation  
  - `Monocle_*.R` — trajectory inference with Monocle3  
  - `SingleR-based_annotation.R` — reference-based cell type annotation  
- `*.csv` — partial processed data tables used for figure reproduction  

---

## Data Availability

Source data: All processed data supporting the findings of this study are provided with the manuscript.

Raw sequencing data and processed h5ad file: Available at the Gene Expression Omnibus (GEO) under accession number GSEXXXXXX.
