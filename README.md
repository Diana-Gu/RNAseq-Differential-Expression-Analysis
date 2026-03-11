# RNA-seq Differential Expression Analysis

This repository presents a basic RNA-seq differential expression workflow implemented in **R** using the **DESeq2** package.
The analysis uses a real RNA-seq dataset (`airway`) and demonstrates common steps used in transcriptomic data analysis.

---

## Project Overview

The goal of this project is to identify genes that are differentially expressed between experimental conditions using RNA-seq count data.

The workflow includes:

* Data loading and preprocessing
* Differential expression analysis with **DESeq2**
* Principal Component Analysis (PCA)
* Heatmap visualization of variable genes
* Volcano plot of differential expression results

---

## Dataset

This analysis uses the **airway RNA-seq dataset**, which contains gene expression data from human airway smooth muscle cells treated with dexamethasone.

The dataset is available through the Bioconductor package:

```
airway
```

---

## Tools and Libraries

The analysis was performed using the following R packages:

* **DESeq2**
* **ggplot2**
* **pheatmap**
* **matrixStats**
* **airway**

---

## Repository Structure

```
RNAseq-Differential-Expression-Analysis
│
├── scripts
│   └── rnaseq_analysis_airway.R
│
├── results
│   ├── deseq2_results.csv
│   ├── pca_plot.png
│   ├── heatmap_top_genes.png
│   └── volcano_plot.png
│
└── README.md
```

---

## How to Run the Analysis

1. Clone this repository:

```
git clone https://github.com/YOUR_USERNAME/RNAseq-Differential-Expression-Analysis.git
```

2. Open **R** or **RStudio**.

3. Install the required packages:

```r
install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("airway")
BiocManager::install("pheatmap")
BiocManager::install("matrixStats")
```

4. Run the analysis script:

```r
source("scripts/rnaseq_analysis_airway.R")
```

The results will be automatically saved in the `results/` directory.

---

## Output

The workflow generates several outputs:

* **Differential expression table** (`deseq2_results.csv`)
* **PCA plot** showing sample clustering
* **Heatmap** of highly variable genes
* **Volcano plot** highlighting significantly regulated genes

---

## Author

**Diana Gutiérrez Martínez**

