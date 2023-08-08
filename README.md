# SpaDo
Cross-slice **Spa**tial Transcriptome **Do**main Analysis.

![](Workflow.png)<!-- -->
**Fig 1. Workflow of SpaDo.** **a.** Calculating the SPACE for both single-cell resolution and spot resolution spatial transcriptomic data. SPACE: SPatially Adjacent Cell type Embedding. **b.** Three functions involved in cross-slice spatial domain analysis: multi-slice domain detection, reference-based spatial domain annotation, and cross-slice clustering analysis by consideration of spatial domain composition. JSD, Jensen Shannon Divergence.
## Overview
With the rapid advancements in spatial transcriptome sequencing, multiple tissue slices become available, allowing for the integration and decoding of spatial cellular landscapes. Herein, we introduce SpaDo, an integrative tool that facilitates cross-slice spatial domain detection, annotation, and downstream analysis for spatial transcriptome at single-cell and spot resolutions. SpaDo contains the following modules: (1) spatial domain detection in a cross-slice way, (2) reference-based spatial domain annotation, and (3) cross-slice clustering analysis. The superiority of SpaDo was demonstrated by a comprehensive investigation on 40 sets of multi-slice spatial transcriptome data from 6 different sequencing platforms. Collectively, our results highlight the tremendous potential of SpaDo to gain novel biological insights across multi-slice spatial transcriptomes.

## Installation
* **SpaDo** package can be installed from Github using **devtools** packages with **R>=4.0.5**.

    ```r
    library(devtools)
    install_github("bm2-lab/SpaDo")
    ```
    
## Getting started
See [Tutorials](https://www.jianguoyun.com/p/DW15NecQnMvoCxji45QFIAA) and [Demo datasets](https://www.jianguoyun.com/p/DX1ssBYQnMvoCxjZ45QFIAA)

## Citation
Bin Duan, Shaoqi Chen, Xiaojie Cheng, Qi Liu. Cross-slice Spatial Transcriptome Domain analysis with SpaDo.

