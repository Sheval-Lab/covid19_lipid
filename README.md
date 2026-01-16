# Lipid Metabolism Dysregulation in the Context of Pneumonia

The repository contains the code to reproduce meta-analysis of bulk and single-cell RNA-seq datasets in the paper "Lung lipid deposition in pneumonias of different etiologies" by Potashnikova D.M. et al. 

## scRNA-seq data reanalysis

- `01_scRNAseq-DGEA.R` - differential gene expression analysis in cell types across 45 scRNA-seq datasets from [SCovid v2.0 database](http://bio-computing.hrbmu.edu.cn/scovid/) ([Zhu et al. 2025](https://doi.org/10.1128/spectrum.01933-24)),
- `02_scRNAseq-ORA.R` over-representation analysis (ORA) based on Gene Ontology and KEGG databases, ORA results keyword-based filtering to look for lipid metabolism-related gene terms and pathways,
- `03_scRNAseq-semicircle-plots.R` - visualizations of enriched lipid metabolism-related gene sets in individual cell types, datasets, or tissue types,
- `04_scRNAseq-log2FC-plots.R` - visualization of average log2 fold change values in individual cell types in datasets of lung and airway tissues.

Directory `results/scRNAseq` contains all results (tables and plots) of this analysis. 