# Spatial Host-Microbiome sequencing (SHM-seq)

Lötstedt B *et al* (2023) Spatial host-microbiome sequencing reveals niches in the mouse gut

Mucosal and barrier tissues such as the gut, lung or skin, are composed of a complex network of cells and microbes forming a tight niche that prevents pathogen colonization and supports host-microbiome symbiosis. Characterizing these networks at high molecular and cellular resolution is crucial for our understanding of homeostasis and disease. Here, we present spatial host-microbiome sequencing (SHM-Seq), an all-sequencing based approach that captures tissue histology, polyadenylated RNAs and bacterial 16S sequences directly from a tissue by modifying spatially barcoded glass surfaces to enable simultaneous capture of host transcripts and hypervariable regions of the 16S bacterial rRNA. We apply our approach to the mouse gut as a model system, use a deep learning approach for data mapping, and detect spatial niches defined by cellular composition and microbial geography. We show that subpopulations of gut cells express specific gene programs in different microenvironments characteristic of regional commensal bacteria and impact host-bacteria interactions. SHM-Seq should enhance the study of native host-microbe interactions in health and disease.

### SHM-seq workflow
![github-small](https://github.com/brittalot/spatial_host_microbiome_sequencing/blob/master/SHM-seq-fig1.png)

Illustration kindly made by Ania Hupalowska.

### Data availability
The raw and processed sequencing and image files needed to recreate all the results in this study have been made avaiable at [Broad's Single Cell Portal]().

### Data pre-processing
Initial host sequncing data processing was performed with ST Pipeline ([v.1.7.6](https://github.com/SpatialTranscriptomicsResearch/st_pipeline/releases/tag/1.7.6)) and initial bacterial sequencing data processing was performed using the [taxonomy assigment pipeline](./7.Bacterial-analysis/bacterial_preparation_scripts/taxonomy_assignment_pipeline/).

For using our spatial spots alignments and reporting tool, please go to our [SpoTteR repository](https://github.com/klarman-cell-observatory/SpoTteR).

### Spatial expression estimates using Splotch
For generating spatial gene expression estimates and spatial differential expression analysis, we advise you to follow instruction at: https://github.com/tare/Splotch and cite Äijö T, Maniatis S & Vickovic S et al: Splotch: Robust estimation of aligned spatial temporal gene expression data, doi: https://doi.org/10.1101/757096. In order to ease use, we have made the complete Splotch workflow available through [Broad's Firecloud platform](https://portal.firecloud.org/?return=firecloud#methods/jgoud/splotch/58).

To use Ilastik for spatial analysis of bacterial fluorescence, please check [pre-processing code](./pre-splotch/) and to implement it with gene expression estimates, please check [post-processing](./post-splotch/).

### Deep learning model used in the taxonomy assignment pipeline
For deep learning model architecture, training and evaluation, please see [DL model](./Bacterial-analysis/bacterial_preparation_scripts/DLmodel/)

  