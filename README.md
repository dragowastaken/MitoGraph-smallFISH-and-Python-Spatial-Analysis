# MitoGraph-smallFISH-and-Python-Spatial-Analysis

<img width="960" height="672" alt="Thesis Figures (5)" src="https://github.com/user-attachments/assets/a4dfb96e-127d-4900-9106-ae0ecb79c19a" />

**Computational workflow for single-cell RNA-mitochondria spatial analysis**

Schematic overview of the image-analysis pipeline used to quantify RNA localization relative to the mitochondrial network. From larger microscopy fields, individual cells were manually selected and cropped so that each cell could be analyzed independently. Cropped cells were separated into individual fluorescence channels for mRNA and mitochondria. RNA puncta coordinates were detected using smallFISH, while mitochondrial morphology was processed using MitoGraph to reconstruct the mitochondrial skeleton. Detected mRNA and mitochondrial skeleton coordinates were converted into a common three-dimensional coordinate system for quantitative analysis. The final output enables comparison of mRNA puncta positions, mitochondrial network morphology, and spatial relationships between RNA species and mitochondria. Scale bars, 5 µm.

<img width="960" height="672" alt="Thesis Figures (6)" src="https://github.com/user-attachments/assets/8e82a7c0-c3a2-4d54-80f6-f6dd81dced00" />

**Representative smFISH image processing and three-dimensional RNA-mitochondria reconstruction.**

Example of the imaging and computational analysis workflow applied to a single S. cerevisiae cell. Representative fluorescence images showing ATP6/8 mRNA, TIM50 mRNA, and mito-GFP-labeled mitochondria from the same cell. RNA signals appear as puncta, while the mitochondrial signal outlines the mitochondrial network. Scale bars, 5 µm. Three-dimensional reconstruction of that same cell showing ATP6/8 mRNA (magenta), TIM50 mRNA (green), and mitochondrial skeleton(grey) in a shared coordinate space. This reconstruction provides the basis for quantitative distance-based analysis of RNA organization relative to mitochondria. 

<img width="960" height="672" alt="Thesis Figures (7)" src="https://github.com/user-attachments/assets/448b270b-ca0c-4361-bbbe-346a7fdb001e" />

**Example WT pooled real and randomized nearest-neighbor distance distributions.**

Pooled nearest-neighbor distance distributions are shown for ATP6/8, ATP2, ATP3, and TIM50 in the WT yWL333 background. Filled histograms represent real detected mRNA nearest-neighbor distances, while blue outline histograms represent matched randomized mRNA distributions. Dashed vertical lines indicate the median distance for each distribution, and Δ median indicates the real-minus-random median shift shown in each plot. Distances below 0.20 µm were floored to 0.20 µm to avoid overinterpreting sub-diffraction-scale distances. 
