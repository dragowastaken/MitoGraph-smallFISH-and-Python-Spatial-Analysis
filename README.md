# MitoGraph-SmallFISH-and-Python-Spatial-Analysis

This repository contains a semi-automated image-analysis workflow developed as part of my thesis work under [Dr. Weihan Li](https://www.weihan-li.com/), submitted in partial fulfillment of the requirements for the **Degree of Master of Science in the Graduate Program of Biotechnology at Brown University**. 

The project was designed to support reproducible, large-scale, single-cell spatial analysis of _Saccharomyces cerevisiae_. It combines [Small Fish](https://github.com/EBL-IGH/small_fish_gui) for smFISH RNA spot detection and quantification with [MitoGraph](https://github.com/vianamp/MitoGraph) for three-dimensional mitochondrial network reconstruction. These outputs are then integrated using custom Python scripts to quantify RNA localization relative to mitochondria across cells, imaging series, probe sets, and experimental conditions.

The workflow can be used to subtract image background when needed, crop raw z-stacks into single-cell images, run MitoGraph and smallFISH, extract RNA intensity distributions, convert mRNA spot coordinates into microns, calculate real mRNA nearest-neighbor distances, test mRNA-mRNA colocalization, measure mRNA-to-mitochondrial-surface proximity, generate randomized mitochondrial-proximity null models, compare real and randomized spatial distributions, and summarize strain/probe effects across experimental conditions.

This repository is intended as both a reproducible analysis pipeline and a practical handoff guide for future users in the Li Lab, as well as others who want to analyze the spatial organization of RNA in mitochondria using microscopy data. A full procedure guide is listed [here](./MitoGraph_smallFISH_Procedure_724).

## Workflow and Important Outputs

<img width="960" height="672" alt="Thesis Figures (5)" src="https://github.com/user-attachments/assets/a4dfb96e-127d-4900-9106-ae0ecb79c19a" />

### **Computational workflow for single-cell RNA-mitochondria spatial analysis**

Schematic overview of the image-analysis pipeline used to quantify RNA localization relative to the mitochondrial network. (A) From larger microscopy fields, individual cells were manually selected and cropped so that each cell could be analyzed independently. (B) Cropped cells were separated into individual fluorescence channels for mRNA and mitochondria. RNA puncta coordinates were detected using SmallFISH, while mitochondrial morphology was processed using MitoGraph to reconstruct the mitochondrial skeleton. (C) Detected mRNA and mitochondrial skeleton coordinates were converted into a common three-dimensional coordinate system for quantitative analysis. The final output enables comparison of mRNA puncta positions, mitochondrial network morphology, and spatial relationships between RNA species and mitochondria. Scale bars, 5 µm.

<img width="960" height="672" alt="Thesis Figures (6)" src="https://github.com/user-attachments/assets/8e82a7c0-c3a2-4d54-80f6-f6dd81dced00" />

### **Representative smFISH image processing and three-dimensional RNA-mitochondria reconstruction.**

Example of the imaging and computational analysis workflow applied to a single S. cerevisiae cell. Representative fluorescence images showing ATP6/8 mRNA, TIM50 mRNA, and mito-GFP-labeled mitochondria from the same cell. RNA signals appear as puncta, while the mitochondrial signal outlines the mitochondrial network. Scale bars, 5 µm. Three-dimensional reconstruction of that same cell showing ATP6/8 mRNA (magenta), TIM50 mRNA (green), and mitochondrial skeleton(grey) in a shared coordinate space. This reconstruction provides the basis for quantitative distance-based analysis of RNA organization relative to mitochondria. Check out an example of [3D Cell Model](./cell_002_GLOBAL_CALIBRATION_surface_overlay.html).

<img width="960" height="672" alt="Thesis Figures (7)" src="https://github.com/user-attachments/assets/448b270b-ca0c-4361-bbbe-346a7fdb001e" />

### **Example WT pooled real and randomized nearest-neighbor distance distributions.**

Pooled nearest-neighbor distance distributions are shown for ATP6/8, ATP2, ATP3, and TIM50 in the WT background. Filled histograms represent real detected mRNA nearest-neighbor distances, while blue outline histograms represent matched randomized mRNA distributions. Dashed vertical lines indicate the median distance for each distribution, and Δ median indicates the real-minus-random median shift shown in each plot. Distances below 0.20 µm were floored to 0.20 µm to avoid overinterpreting sub-diffraction-scale distances. 

## Limitations and Future Work: 
- Image quality must be high for reliable downstream analysis. Additional image-improvement methods, such as deconvolution, have not yet been explored extensively.
- Selection of high-quality single cells currently requires substantial manual input, especially when analyzing statistically meaningful sample sizes of more than 100 cells. A useful future extension would be automated single-cell detection, filtering, and cropping.

## [Example Data Set ](https://drive.google.com/drive/folders/1GKXuFGcioTF_7YcUi1d277WfLsQV1za7?usp=sharing)

## Acknowledgments

This project was completed as part of my thesis work in the Li Lab at Brown University. I am especially grateful to my advisor, [Dr. Weihan Li](https://www.weihan-li.com/), for his guidance, mentorship, and support throughout the development of this project.

I am particularly grateful to Dubuke Ma, whose foundational code helped support the computational analysis pipeline.

This workflow also benefited from AI-assisted coding and documentation support. [ChatGPT Codex](https://learn.chatgpt.com/docs) was used as a development assistant to help accelerate script writing, debugging, documentation, and organization of the analysis pipeline. All analysis goals, biological interpretation, workflow decisions, and final implementation choices were directed and reviewed by the author.

This work was supported by the NIH Common Fund, grant R00GM148788. Figures included were created with [BioRender](https://www.biorender.com/).

<img width="82" height="20" alt="image" src="https://github.com/user-attachments/assets/11749e98-2c49-4788-bfda-055d39ee4c28" />

