# DLBCL-TPS12

## Single-cell dissection reveals a conserved stress-adaptive proliferative program in diffuse large B-cell lymphoma

Repository containing all analysis scripts used for the study of TPS12, a biologically derived transcriptional program identified from malignant-cell states in diffuse large B-cell lymphoma.

------------------------------------------------------------------------

## Overview

This repository contains the R scripts used to identify, derive, and validate TPS12, a biologically derived transcriptional program associated with stress-adaptive proliferative malignant-cell states in diffuse large B-cell lymphoma (DLBCL).

The study integrates single-cell and bulk transcriptomic datasets to:

-   derive TPS12 from malignant-cell states
-   validate TPS12 across independent cohorts
-   evaluate survival associations
-   characterize longitudinal remodeling during relapse
-   compare TPS12 with a published 30-gene relapse discriminator

## Graphical abstract

![Graphical Abstract](Graphical%20abstract.png)

---

## Public datasets

| Dataset | Purpose |
|---------|---------|
| GSE182434 | Single-cell discovery |
| GSE31312 | Bulk and molecular-context validation |
| GSE10846 | Survival validation |
| GSE32918 | External validation |
| GSE23501 | Independent COO validation |
| GSE193566 | Longitudinal diagnosis-relapse analysis |
| NCI DLBCL cohort | Validation across published genetic subtypes |
## Workflow

1. Single-cell preprocessing
2. Malignant-cell identification
3. TPS12 derivation
4. Bulk validation
5. Survival analysis
6. Molecular-context validation across COO and published genetic subtypes
7. Longitudinal validation
8. Comparison with published relapse discriminator
## Software

R 4.6

Main packages

-   Seurat
-   CopyKAT
-   limma
-   GSVA
-   survival
-   survminer
-   tidyverse
-   ggplot2
-   patchwork

## Repository structure

scripts/ figures/ results/ data/ README.md LICENSE

## Citation

If you use this repository, please cite the accompanying manuscript.

Ahmed Galal Rezk

Single-cell dissection reveals a conserved stress-adaptive proliferative program in diffuse large B-cell lymphoma.


