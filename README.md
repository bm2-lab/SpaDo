# [SpaDo](https://doi.org/10.1186/s13059-024-03213-x)
Multi-slice **Spa**tial Transcriptome **Do**main Analysis.

![](Overview.png)<!-- -->
**Fig 1. Workflow of SpaDo.** **a.** Calculating the SPACE for both single-cell resolution and spot resolution spatial transcriptomic data. SPACE: SPatially Adjacent Cell type Embedding. **b.** Three functions involved in multi-slice spatial domain analysis: multi-slice domain detection, reference-based spatial domain annotation, and multi-slice clustering analysis by consideration of spatial domain composition. JSD, Jensen Shannon Divergence.
## Overview
With the rapid advancements in spatial transcriptome sequencing, multiple tissue slices are now available, enabling the integration and interpretation of spatial cellular landscapes. Herein, we introduce SpaDo, a tool for multi-slice spatial domain analysis, including modules for multi-slice spatial domain detection, reference-based annotation, and multiple slice clustering at both single-cell and spot resolutions. We demonstrate SpaDo's effectiveness with over 40 multi-slice spatial transcriptome datasets from 7 sequencing platforms. Our findings highlight SpaDo's potential to reveal novel biological insights in multi-slice spatial transcriptomes.

## Update (2024.8.9)
* **(1)** SpaDo was originally built on **Seurat 4.0.4** and **SeuratObject 4.0.4**. To avoid problems due to Seurat updates, We recently (2024.8.9) updated SpaDo to ensure compatibility with **Seurat v5**;
* **(2)** Enabling parallel execution for the **"DistributionDistance"** function when using the "JSD" distance metric, significantly reducing runtime;
* **(3)** We updated corresponding [Tutorials](https://www.jianguoyun.com/p/DW15NecQnMvoCxji45QFIAA).

## Installation
* **SpaDo** package can be installed from Github using **devtools** packages with **R>=4.0.5**.

    ```r
    library(devtools)
    install_github("bm2-lab/SpaDo")
    ```
   
## Getting started
Here are [Tutorials](https://www.jianguoyun.com/p/DW15NecQnMvoCxji45QFIAA) with different styles and [Demo datasets](https://www.jianguoyun.com/p/DX1ssBYQnMvoCxjZ45QFIAA)
used in our tutorials.
## Citation
Duan, B., Chen, S., Cheng, X. et al. Multi-slice spatial transcriptome domain analysis with SpaDo. Genome Biol 25, 73 (2024). https://doi.org/10.1186/s13059-024-03213-x

