# MitoGraph, smallFISH, and Python Spatial-Analysis Procedure

**Purpose of this document:** This is a folder-level handoff guide for the full Li Lab image-analysis workflow. A new user should be able to open an experiment folder, understand what was run, identify the important inputs and outputs, and trace how each script feeds the next step.

The workflow subtracts image background when needed, converts raw z-stacks into single-cell images, runs MitoGraph and smallFISH, extracts raw smallFISH RNA intensity distributions, converts mRNA spot coordinates into microns, computes real mRNA nearest-neighbor distributions, checks mRNA-mRNA colocalization cell-by-cell, measures optional mRNA-to-mitochondrial-surface proximity, builds random mitochondrial-proximity null models, compares real versus random distributions, and finally summarizes strain/probe effects across experimental conditions.

---

## Table of Contents

This table of contents is intentionally placed at the beginning so someone opening the folder can quickly jump to the relevant workflow stage or output family.

- [1. Core concept](#1-core-concept)
- [2. Workflow summary](#2-workflow-summary)
- [3. Calibration and coordinate assumptions](#3-calibration-and-coordinate-assumptions)
- [4. Expected folder layout after image preprocessing](#4-expected-folder-layout-after-image-preprocessing)
- [PART A — ImageJ/Fiji preprocessing macros](#part-a-imagejfiji-preprocessing-macros)
  - [4A. `0.Recursive_Subtract_Background_OVERWRITE_ORIGINALS.ijm`](#4a-0recursive_subtract_background_overwrite_originals1ijm)
  - [5. `1.Generate_multi_channel_MaxProjs_per_subfolder_GRAYSCALE_AUTOSCALE.ijm`](#5-1generate_multi_channel_maxprojs_per_subfolder_grayscale_autoscaleijm)
  - [6. `2.CropCells_recursive_per_subfolder_MATCHFIX_v2.ijm`](#6-2cropcells_recursive_per_subfolder_matchfix_v2ijm)
  - [7. `3.Stack to Hyperstack_recursive_cells_v2.ijm`](#7-3stack-to-hyperstack_recursive_cells_v2ijm)
  - [8. `4.Extract Channels_recursive_Hyperstacks_Grayscale_named_channels.ijm`](#8-4extract-channels_recursive_hyperstacks_grayscale_named_channelsijm)
- [PART B — MitoGraph](#part-b-mitograph)
  - [9. Running MitoGraph](#9-running-mitograph)
- [PART C — smallFISH](#part-c-smallfish)
  - [10. Running smallFISH](#10-running-smallfish)
  - [10A. smallFISH RNA intensity histogram QC](#10a-smallfish-rna-intensity-histogram-qc)
- [PART D — Python scripts](#part-d-python-scripts)
  - [11. `5.mito_vis_per_cells_folder_recursive.py`](#11-5mito_vis_per_cells_folder_recursivepy)
  - [12. `6.mrna_nn_distance_recursive.py`](#12-6mrna_nn_distance_recursivepy)
  - [12A. `rna_intensity_histograms_recursive.py`](#12a-rna_intensity_histograms_recursivepy)
  - [13. `mrna_colocalization_by_cell.py`](#13-mrna_colocalization_by_cellpy)
  - [14. `7.random_mrna_null_recursive.py`](#14-7random_mrna_null_recursivepy)
  - [15. `7.1 random_mrna_distance_range.py`](#15-71-random_mrna_distance_rangepy)
  - [16. `8.real_vs_random_recursive_comparison.py`](#16-8real_vs_random_recursive_comparisonpy)
  - [17. `9.pooled_real_vs_random_across_series_REPLICATE_LEVEL_STATS.py`](#17-9pooled_real_vs_random_across_series_replicate_level_statspy)
  - [18. `10. pooled_real_vs_random_ Thesis Graphs.py`](#18-10-pooled_real_vs_random_-thesis-graphspy)
  - [19. `11.compare_strains_probe_sets_boxplot_v6.py`](#19-11compare_strains_probe_sets_boxplot_v6py)
- [PART E — Optional mitochondrial proximity and node-distance scripts](#part-e-optional-mitochondrial-proximity-and-node-distance-scripts)
  - [20. `mito_rna_surface_GLOBAL_CALIBRATION.py`](#20-mito_rna_surface_global_calibrationpy)
  - [20A. `mito_rna_surface_CORRECTED_XY_0p065_CROP250.py`](#20a-mito_rna_surface_corrected_xy_0p065_crop250py)
  - [21. `mrna_to_nearest_node_distance_BY_SERIES_FIXED.py`](#21-mrna_to_nearest_node_distance_by_series_fixedpy)
  - [22. `mrna_to_nearest_node_distance_analysis.py` and older `mrna_to_nearest_node_distance_BY_SERIES.py`](#22-mrna_to_nearest_node_distance_analysispy-and-older-mrna_to_nearest_node_distance_by_seriespy)
- [PART F — Output map: what to look for in a completed folder](#part-f-output-map-what-to-look-for-in-a-completed-folder)
  - [23. At the Series level](#23-at-the-series-level)
  - [24. At the condition parent level](#24-at-the-condition-parent-level)
  - [25. At the experiment parent level across strains/probes](#25-at-the-experiment-parent-level-across-strainsprobes)
- [PART G — How the scripts connect](#part-g-how-the-scripts-connect)
  - [26. Dependency chain](#26-dependency-chain)
  - [27. Minimum files required for a completed real-vs-random condition](#27-minimum-files-required-for-a-completed-real-vs-random-condition)
- [PART H — QC checklist](#part-h-qc-checklist)
  - [27A. After background subtraction](#27a-after-background-subtraction)
  - [28. Before MitoGraph](#28-before-mitograph)
  - [29. After MitoGraph](#29-after-mitograph)
  - [30. After smallFISH](#30-after-smallfish)
  - [30A. After RNA intensity histogram script](#30a-after-rna-intensity-histogram-script)
  - [31. After real NN script](#31-after-real-nn-script)
  - [32. After mRNA-mRNA colocalization script](#32-after-mrna-mrna-colocalization-script)
  - [33. After mitochondrial surface calibration script](#33-after-mitochondrial-surface-calibration-script)
  - [34. After random null script](#34-after-random-null-script)
  - [35. After pooled analysis](#35-after-pooled-analysis)
  - [36. After cross-strain comparison](#36-after-cross-strain-comparison)
- [PART I — Common failure modes and fixes](#part-i-common-failure-modes-and-fixes)
  - [37. No `spots_extraction` folders found](#37-no-spots_extraction-folders-found)
  - [37A. RNA intensity histogram script finds no intensity values](#37a-rna-intensity-histogram-script-finds-no-intensity-values)
  - [37B. Background subtraction is too weak, too strong, or accidentally repeated](#37b-background-subtraction-is-too-weak-too-strong-or-accidentally-repeated)
  - [38. Real NN script makes no NN distances](#38-real-nn-script-makes-no-nn-distances)
  - [39. Random null script reports missing shared skeletons](#39-random-null-script-reports-missing-shared-skeletons)
  - [40. Real-vs-random script cannot find random files](#40-real-vs-random-script-cannot-find-random-files)
  - [41. mRNA colocalization script finds no matched cells](#41-mrna-colocalization-script-finds-no-matched-cells)
  - [42. Surface calibration script reports missing RNA files or missing mito surfaces](#42-surface-calibration-script-reports-missing-rna-files-or-missing-mito-surfaces)
  - [42A. Corrected surface overlay is mirrored or shifted](#42a-corrected-surface-overlay-is-mirrored-or-shifted)
  - [43. Windows path-length errors](#43-windows-path-length-errors)
  - [44. Node-distance script matches only `cell_002` or mismatches `_000_nodes.vtk`](#44-node-distance-script-matches-only-cell_002-or-mismatches-_000_nodesvtk)
  - [45. Derived surface-distance CSVs are treated as input mRNA files](#45-derived-surface-distance-csvs-are-treated-as-input-mrna-files)
  - [46. Using different probes other than MS2, ATP2, ATP3, or TIM50](#46-using-different-probes-other-than-ms2-atp2-atp3-or-tim50)
- [PART J — Final handoff notes](#part-j-final-handoff-notes)

---

## 1. Core concept

The analysis is organized around **one parent folder per experimental condition**. An experimental condition is usually a strain/probe/treatment folder containing multiple imaging series.

Typical condition folder:

```text
Condition parent folder/
├── Series 4/
│   ├── cells/
│   ├── MS2/
│   ├── ATP2/ or ATP3/ or TIM50/
│   └── ...
├── Series 5/
├── Series 8/
├── Series 9/
└── Series 11/
```

Most Python scripts ask you to select this parent folder with a folder-picker dialog. This avoids hard-coded paths and makes the workflow reusable across strains and probe sets.

The downstream logic assumes that files are matched **within the same Series folder**. Do not mix mRNA files from one series with MitoGraph skeletons or node files from another series.

---

## 2. Workflow summary

Run the workflow in this order:

0. `0.Recursive_Subtract_Background_OVERWRITE_ORIGINALS.ijm`  
   Optional but recommended raw-image cleanup step before projection/cropping. It recursively subtracts background from TIFF stacks and overwrites the originals. The standard rolling-ball radius used in this workflow is `30` pixels. Make a backup before running.

1. `1.Generate_multi_channel_MaxProjs_per_subfolder_GRAYSCALE_AUTOSCALE.ijm`  
   Creates quick multi-channel maximum projections for ROI selection. Output: `MaxProjs.tif`.

2. `2.CropCells_recursive_per_subfolder_MATCHFIX_v2.ijm`  
   Crops individual cells using `RoiSet.zip`. Output: `cells/*.tif`.

3. `3.Stack to Hyperstack_recursive_cells_v2.ijm`  
   Converts cropped stacks into proper hyperstacks. Output: `cells/Hyperstacks_Grayscale/*.tif`.

4. `4.Extract Channels_recursive_Hyperstacks_Grayscale_named_channels.ijm`  
   Splits selected channels into named channel folders. Output: `Hyperstacks_Grayscale/<ChannelName>/*.tif`.

5. MitoGraph  
   Skeletonizes mitochondrial channel images. Output: MitoGraph `.txt`, `_nodes.vtk`, `_mitosurface.vtk`, and related files.

6. smallFISH  
   Detects mRNA spots. Output: `results/spots_extraction/spots_extractions_*.csv` or `.xlsx`.

6A. `rna_intensity_histograms_recursive.py`  
   Optional smallFISH intensity-QC script. Extracts the raw `intensity` column from smallFISH spot tables and plots RNA intensity histograms by Series, experimental condition, and RNA channel. Output folder: `RNA_intensity_histograms/`.

7. `5.mito_vis_per_cells_folder_recursive.py`  
   QC of MitoGraph skeletons. Output: `cells/mito_visualization_interactive/*.html`.

8. `6.mrna_nn_distance_recursive.py`  
   Converts smallFISH coordinates and computes real mRNA nearest-neighbor distances. Outputs: `converted_coordinates/`, `*_NN_distance_um.npy`, `*_NN_distance_um.xlsx`, `*_spot_counts.csv`, and `*_spot_counts.xlsx`.

9. `mrna_colocalization_by_cell.py`  
   Optional mRNA-mRNA colocalization analysis. Matches two mRNA channels cell-by-cell within each Series and reports puncta within an adjustable 3D distance cutoff. Output folder: `mrna_coloc_BY_SERIES/`.

10. `mito_rna_surface_GLOBAL_CALIBRATION.py` or `mito_rna_surface_CORRECTED_XY_0p065_CROP250.py`  
    Optional mRNA-to-mitochondrial-surface proximity analysis. Use the corrected `CROP250` version for datasets processed with the corrected MitoGraph microscope calibration. Output folder: `mito_rna_surface_CORRECTED_XY_0p065_CROP250/`.

11. `7.random_mrna_null_recursive.py`  
    Generates the standard random mRNA null model near MitoGraph skeletons. Output: `random_<mRNA>_output/random_<mRNA>_nn.npy`.

12. `7.1 random_mrna_distance_range.py`  
    Optional sensitivity/null model with a custom skeleton-distance range. Output: distance-tagged `random_<mRNA>_output_dist_*/*_dist_*.npy` files.

13. `8.real_vs_random_recursive_comparison.py`  
    Compares real vs random within each series/folder. Output: `plots/<mRNA>_real_vs_random.png`.

14. `9.pooled_real_vs_random_across_series_REPLICATE_LEVEL_STATS.py`  
    Pools all series for one condition with replicate-level statistics. Output folder: `pooled_real_vs_random_comparisons/`.

15. `10. pooled_real_vs_random_ Thesis Graphs.py`  
    Thesis-figure version of pooled plots with stricter visual inference gates. Output: strain-labeled pooled plots and CSVs.

16. `11.compare_strains_probe_sets_boxplot_v6.py`  
    Compares strains and probe sets across condition folders. Output folder: `strain_probe_comparison/`.

17. `mrna_to_nearest_node_distance_BY_SERIES_FIXED.py`  
    Optional mRNA-to-MitoGraph-node distance analysis. Output folder: `node_mrna_dist_BY_SERIES/`.

## 3. Calibration and coordinate assumptions

### MitoGraph calibration

MitoGraph should be run with the microscope voxel calibration used for the mitochondrial channel. The common command used in this workflow is:

```powershell
.\MitoGraph.exe -xy 0.0645 -z 0.2 -path "C:\path\to\Series XX\cells"
```

Some scripts use the more exact XY value:

```python
X_SCALE = 0.0645
Y_SCALE = 0.0645
Z_SCALE = 0.2
```

Before starting a new dataset, *verify the microscope calibration and keep it consistent across MitoGraph, smallFISH coordinate conversion, and downstream Python scripts.*

### smallFISH coordinate order

The real mRNA NN script assumes smallFISH spot coordinates are stored as:

```text
[z, y, x] in pixels
```

It converts them into:

```text
[x_um, y_um, z_um]
```

using the scale factors above.

---

## 4. Expected folder layout after image preprocessing

A completed condition folder usually contains:

```text
Condition parent folder/
├── Series 4/
│   ├── MaxProjs.tif
│   ├── RoiSet.zip
│   ├── cells/
│   │   ├── cropped single-cell TIFFs
│   │   ├── MitoGraph output .txt files
│   │   ├── MitoGraph *_nodes.vtk files
│   │   ├── Hyperstacks_Grayscale/
│   │   └── mito_visualization_interactive/
│   ├── MS2/
│   │   └── .../results/spots_extraction/
│   │       ├── spots_extractions_*.csv or .xlsx
│   │       ├── converted_coordinates/
│   │       ├── MS2_NN_distance_um.npy
│   │       ├── MS2_NN_distance_um.xlsx
│   │       ├── MS2_spot_counts.csv
│   │       ├── MS2_spot_counts.xlsx
│   │       ├── random_MS2_output/
│   │       └── plots/
│   └── ATP2/ ATP3/ or TIM50/
│       └── .../results/spots_extraction/
├── RNA_intensity_histograms/               optional, if RNA intensity QC was run
├── pooled_real_vs_random_comparisons/
├── mrna_coloc_BY_SERIES/                    optional, if mRNA-mRNA colocalization was run
├── mito_rna_surface_GLOBAL_CALIBRATION/     older optional surface-proximity output
├── mito_rna_surface_CORRECTED_XY_0p065_CROP250/  corrected optional surface-proximity output
└── node_mrna_dist_BY_SERIES/                optional, if node-distance analysis was run
```

---

# PART A — ImageJ/Fiji preprocessing macros


## 4A. `0.Recursive_Subtract_Background_OVERWRITE_ORIGINALS.ijm`

### Purpose

Recursively subtracts image background from TIFF stacks before the rest of the Fiji/ImageJ preprocessing workflow. This helps reduce diffuse fluorescence background that can add noise to both smallFISH spot detection and MitoGraph mitochondrial segmentation.

This is usually run before making maximum projections, drawing ROIs, cropping cells, extracting channels, running smallFISH, or running MitoGraph.

### Important warning

This macro **overwrites the original TIFF files**. Make a backup copy of the raw data before running it. Do not run it repeatedly on the same images unless that is intentional, because each run subtracts background again.

### Input

```text
Any parent folder containing .tif or .tiff image stacks
```

The macro searches recursively through all subfolders under the selected parent folder.

### Standard setting used in this workflow

Use the default rolling-ball radius:

```text
Rolling ball radius = 30 pixels
```

This value is used as the standard because it removes broad background signal without intentionally removing the puncta/skeleton signal that the downstream programs need.

### Dialog settings

The macro opens a settings dialog with these defaults:

```text
Rolling ball radius (pixels): 30
Light background: false
Separate colors: false
Create background only: false
Sliding paraboloid: false
Disable smoothing: false
Process stack slices: true
```

For the standard workflow, keep the default `30` radius and `Process stack slices = true`.

### Output

There is no separate output folder. The macro saves the processed TIFF to the same path, overwriting the original file.

Example:

```text
Before:
Series 4/raw_image.tif

After:
Series 4/raw_image.tif    same filename, background-subtracted image
```

### Why this matters downstream

Diffuse background can cause smallFISH and MitoGraph to detect structures that are not true RNA puncta or mitochondrial signals. Background subtraction reduces this background contribution before:

- `MaxProjs.tif` creation,
- ROI cropping,
- channel extraction,
- smallFISH spot detection,
- MitoGraph skeleton/surface generation.

### QC recommendation

After running this macro, open a few representative images in Fiji and compare them to the backup raw images. Check that:

- real puncta are still visible,
- mitochondria are still continuous enough for MitoGraph,
- background is reduced but not over-flattened,
- no image looks black/clipped,
- the macro was not accidentally run twice on the same dataset.

---

## 5. `1.Generate_multi_channel_MaxProjs_per_subfolder_GRAYSCALE_AUTOSCALE.ijm`

### Purpose

Creates one grayscale, autoscaled, multi-channel maximum-intensity projection stack in each folder that directly contains TIFF images.

### Input

```text
Any folder under the selected parent that directly contains .tif or .tiff z-stacks
```

The macro skips files already named `MaxProjs.tif` or `MaxProjs.tiff`.

### Output

```text
<folder containing TIFFs>/MaxProjs.tif
```

### What the output is used for

`MaxProjs.tif` is used for visual inspection and ROI drawing before cropping individual cells. 
To draw ROI's, use the ROI manager in Fiji with the rectangle tool to select single cells and add them to the ROI manager, then save it to the parent folder. It should be saved alongside `RoiSet.zip` before running the crop macro.

---

## 6. `2.CropCells_recursive_per_subfolder_MATCHFIX_v2.ijm`

### Purpose

Recursively searches for folders containing both:

```text
MaxProjs.tif
RoiSet.zip
```

It then crops each ROI from the matching raw z-stack image to create individual single-cell TIFF stacks. *All downstream analysis assumes that the pixels are set to 200x200*

### Input

```text
Series folder/
├── MaxProjs.tif
├── RoiSet.zip
└── original TIFF z-stack files
```

### Output

```text
Series folder/cells/<original_file>_<cell_index>.tif
```

Example:

```text
Series 4/cells/MS2 + ATP3..._000.tif
Series 4/cells/MS2 + ATP3..._001.tif
```

### Notes

The macro was adjusted to improve filename matching between ROI-derived names and source TIFFs. The `cells` folder is the main input folder for the MitoGraph and hyperstack/channel extraction steps.

---

## 7. `3.Stack to Hyperstack_recursive_cells_v2.ijm`

### Purpose

Recursively finds folders named `cells` and converts cropped TIFF stacks into grayscale hyperstacks with user-specified channel, z-slice, frame, and stack-order settings.

### Input

```text
Series XX/cells/*.tif or *.tiff
```

### Output

```text
Series XX/cells/Hyperstacks_Grayscale/<same_filename>.tif
```

### What the output is used for

The output hyperstacks are used by the channel-extraction macro, which separates the mitochondrial channel from the smFISH channels.

---

## 8. `4.Extract Channels_recursive_Hyperstacks_Grayscale_named_channels.ijm`

### Purpose

Recursively searches for `Hyperstacks_Grayscale` folders, splits hyperstack channels, and saves selected channels into user-named folders.

### Input

```text
Series XX/cells/Hyperstacks_Grayscale/*.tif
```

### User choices

You enter:

- channel number or numbers to extract,
- a descriptive output folder name for each selected channel.

Example channel folder names:

```text
cells (mitoGFP channel has to be named as cells as folder name for MitoGraph.exe to recognize)
MS2
ATP2
ATP3
TIM50
DIC
```

### Output

```text
Series XX/cells/Hyperstacks_Grayscale/<ChannelName>/<original_filename>.tif
```

### What the output is used for

- Mitochondrial channel images are used for MitoGraph extracting mito skeleton (HAVE TO BE NAMED cells).
- mRNA channel images are used for smallFISH.

### Notes
**Move the extracted channel folders back to the original series folder**; if done correctly, it should ask to replace X files; click yes to replace (this gets rid of the useless stacks produced by MitoGraph). The final file structure should be: 

Condition parent folder/
├── Series 4/
│   ├── MaxProjs.tif
│   ├── RoiSet.zip
│   ├── cells (The Mito channel only)
│   │    ├── Hyperstacks_Grayscale/
│   ├── MS2
│   └── ATP2/ ATP3/ or TIM50/

---

# PART B — MitoGraph

## 9. Running MitoGraph

### Preparation
Go to parent folder --> series X --> cells --> Hyperstacks_Grayscale --> Extracted Channels 
Select the "cells" "MS2" "ATP2" folders and cut and paste them into the "series X" layer to overwrite the previous cropped, unsplit images
Now, the "cells" folder is ready for mitograph recognition

### Purpose

MitoGraph reconstructs the mitochondrial network from single-cell mitochondrial images. Its skeleton coordinate files are used for random mRNA null-model generation and for optional mRNA-to-node analyses.

### Basic command

```powershell
cd "C:\path\to\MitoGraph"
.\MitoGraph.exe -xy 0.0645 -z 0.2 -path "C:\path\to\Series XX\cells"
```

### Recursive PowerShell pattern

Use this when one parent folder contains multiple `Series XX/cells` folders:

```powershell
$parent = "C:\path\to\condition_parent"
$mitograph = "D:\Michael Thesis Image Organization\Scripts\MitoGraph.exe"
Get-ChildItem -Path $parent -Directory -Recurse |
Where-Object { $_.Name -eq "cells" } |
ForEach-Object {
    Write-Host "Running MitoGraph on $($_.FullName)"
    & $mitograph -xy 0.0645 -z 0.2 -path "$($_.FullName)"
}

```

### Expected outputs

MitoGraph writes outputs into the same `cells` folder. The important downstream files are:

```text
Series XX/cells/*.txt
Series XX/cells/*_nodes.vtk
```

The `.txt` skeleton files should contain at least:

```text
x
y
z
width_(um)
```

The `*_nodes.vtk` files are used by the optional node-distance scripts.
The user can also visualize the mito skeleton and mito surface by dragging the nodes.vtk skeleton.vtk and mitosurface.vtk files into Paraview.exe

### QC recommendation

Run the MitoGraph visualization script before using MitoGraph outputs for randomization. *Visually check that the skeletons are centered on mitochondrial structures and do not contain obvious segmentation failures.*

---

# PART C — smallFISH

smallFISH installation guide: https://pypi.org/project/small-fish-gui/

## 10. Running smallFISH

### Purpose

smallFISH detects mRNA spots from smFISH channel images. These spot coordinate tables become the input for the real mRNA nearest-neighbor script.

### Launch

```bash
conda activate small_fish
python -m small_fish_gui
```

### Batch detection settings used in this workflow

| Parameter | Setting |
|---|---|
| Coordinates | `Z=0, Y=1, X=2` |
| Threshold | usually around `50–70`; commonly `50` |
| Voxel size | `1, 2, 3` |
| Spot size | `2, 3, 4` |
| Name Batch| MS2/ATP2/ATP3 etc.. |
### Expected output

Each mRNA channel folder should contain a smallFISH results folder similar to:

```text
Series XX/MS2/MS2_<date-time>/results/spots_extraction/
Series XX/ATP2/ATP2_<date-time>/results/spots_extraction/
Series XX/ATP3/ATP3_<date-time>/results/spots_extraction/
Series XX/TIM50/TIM50_<date-time>/results/spots_extraction/
```

Inside each `spots_extraction` folder, smallFISH writes cell-level coordinate tables, for example:

```text
spots_extractions_<image_name>_000.csv
spots_extractions_<image_name>_001.csv
...
```

These are the direct inputs for `6.mrna_nn_distance_recursive.py` and can also be used for RNA-intensity QC before coordinate conversion.

## 10A. smallFISH RNA intensity histogram QC

### Purpose

After smallFISH finishes spot detection, the raw `spots_extractions_*.csv` or `.xlsx` files contain an `intensity` column for each detected RNA punctum. The optional RNA-intensity histogram script uses this column to check whether the raw spot-intensity distributions look reasonable across Series, conditions, and RNA channels.

This is a **smallFISH QC and comparison step**. It does not alter the raw smallFISH files, does not convert coordinates, and does not feed directly into the real-vs-random nearest-neighbor pipeline.

Use this step to ask questions such as:

```text
Are RNA spot intensities similar across Series within the same condition?
Does one experimental condition have systematically brighter or dimmer detected RNA spots?
Are MS2 and ATP2/ATP3/TIM50 intensity distributions comparable enough for downstream interpretation?
```

### Input

The script reads raw smallFISH output files:

```text
Condition parent folder/
└── Series XX/
    └── <mRNA>/.../results/spots_extraction/
        ├── spots_extractions_..._000.csv or .xlsx
        ├── spots_extractions_..._001.csv or .xlsx
        └── ...
```

The script is designed to handle the semicolon-delimited CSV format produced by smallFISH, for example files with columns such as:

```text
spots_id;spot_id;axis-0;axis-1;axis-2;intensity;coordinates
```

The important column is:

```text
intensity
```

### Script

Use:

```text
rna_intensity_histograms_recursive.py
```

Select the parent folder that contains the Series folders or the larger experiment folder containing multiple condition folders.

### Key settings

```python
HISTOGRAM_BINS = 50
INTENSITY_XMIN = None
INTENSITY_XMAX = None
MAKE_LOG10_HISTOGRAMS = True
FILTER_MIN_INTENSITY = None
FILTER_MAX_INTENSITY = None
```

If `INTENSITY_XMIN` and `INTENSITY_XMAX` are left as `None`, each histogram uses automatic x-axis limits. To force all plots to use the same x-axis, set values such as:

```python
INTENSITY_XMIN = 0
INTENSITY_XMAX = 3000
```

### Output folder

```text
Condition or experiment parent folder/RNA_intensity_histograms/
```

### Main outputs

```text
RNA_intensity_histograms/
├── ALL_RNA_spot_intensities.csv
├── ALL_RNA_spot_intensities.xlsx
├── RNA_intensity_histogram_analysis_workbook.xlsx
├── RNA_intensity_file_summary.csv
├── RNA_intensity_file_summary.xlsx
├── RNA_intensity_summary_by_series.csv
├── RNA_intensity_summary_by_series.xlsx
├── RNA_intensity_summary_by_series_and_channel.csv
├── RNA_intensity_summary_by_series_and_channel.xlsx
├── RNA_intensity_summary_by_condition.csv
├── RNA_intensity_summary_by_condition.xlsx
├── RNA_intensity_summary_by_condition_and_channel.csv
├── RNA_intensity_summary_by_condition_and_channel.xlsx
├── RNA_intensity_summary_by_channel.csv
├── RNA_intensity_summary_by_channel.xlsx
├── RNA_intensity_plot_index.csv
├── RNA_intensity_plot_index.xlsx
├── RNA_intensity_rejected_files.csv
├── RNA_intensity_rejected_files.xlsx
├── RNA_intensity_RUN_SETTINGS.csv
└── RNA_intensity_RUN_SETTINGS.xlsx
```

### Plot outputs

```text
RNA_intensity_histograms/plots/
├── parent/
│   ├── ALL_PARENT_RNA_intensity_histogram.png
│   └── ALL_PARENT_RNA_intensity_histogram_log10.png
├── by_series/
├── by_series_and_channel/
├── by_condition/
└── by_condition_and_channel/
```

Each plot folder contains regular intensity histograms and, when enabled, log10 intensity histograms.

### Output interpretation

- `ALL_RNA_spot_intensities.xlsx`: one row per detected RNA spot, including condition label, Series label, RNA channel, cell index, source file, and intensity value.
- `RNA_intensity_summary_by_series.xlsx`: Series-level summary statistics pooled across RNA channels.
- `RNA_intensity_summary_by_series_and_channel.xlsx`: Series-level summary statistics split by RNA channel.
- `RNA_intensity_summary_by_condition.xlsx`: condition-level summary statistics pooled across RNA channels.
- `RNA_intensity_summary_by_condition_and_channel.xlsx`: condition-level summary statistics split by RNA channel.
- `RNA_intensity_file_summary.xlsx`: audit table showing how many intensity values were extracted from each raw smallFISH file.
- `RNA_intensity_plot_index.xlsx`: index of all histogram PNG files created.
- `RNA_intensity_rejected_files.xlsx`: files that were found but could not be processed or were rejected.
- `RNA_intensity_RUN_SETTINGS.xlsx`: exact settings used for reproducibility.

### Recommended use

Use the intensity histograms as a QC layer before interpreting downstream spatial results. Large shifts in intensity distribution between Series or conditions may indicate differences in image acquisition, thresholding, spot-detection behavior, probe performance, or sample quality. These intensity plots should not replace spatial analyses, but they help explain whether downstream differences could be influenced by detection-intensity differences.

---


# PART D — Python scripts

All Python scripts use a folder-picker dialog near the top of the script. Select the parent folder requested by the script. In most cases this is the condition folder that contains all `Series XX` folders.

Recommended packages:

```bash
C:\ProgramData\miniconda3\python.exe -m pip install numpy pandas scipy matplotlib plotly openpyxl
```

---

## 11. `5.mito_vis_per_cells_folder_recursive.py`

### Purpose

Creates interactive 3D HTML plots of MitoGraph skeleton `.txt` files.

### Input

```text
Condition parent folder/
└── Series XX/
    └── cells/
        └── MitoGraph .txt files with x, y, z, width_(um)
```

### Output

Inside each `cells` folder:

```text
Series XX/cells/mito_visualization_interactive/*.html
```

### Output meaning

Each HTML file is an interactive 3D scatter plot. Marker size and color are based on mitochondrial width. Use these files for MitoGraph QC before accepting the skeletons for randomization.

### Run summary printed to console

The script prints:

```text
Cells folders processed
HTML files created
Files skipped/errors
```

No parent-level CSV is written by this visualization script.

---

## 12. `6.mrna_nn_distance_recursive.py`

### Purpose

Converts smallFISH coordinate tables into microns and calculates **real mRNA nearest-neighbor distances** within each mRNA channel.

For every detected mRNA spot in a cell, the script finds the nearest other spot of the same mRNA in that same cell. These real mRNA-to-mRNA nearest-neighbor distributions are used later for real-vs-random comparisons.

### Input

```text
Condition parent folder/
└── Series XX/
    └── <mRNA>/.../results/spots_extraction/
        ├── spots_extractions_..._000.csv or .xlsx
        ├── spots_extractions_..._001.csv or .xlsx
        └── ...
```

Supported mRNA names are inferred from paths and filenames, including:

```text
MS2, ATP2, ATP3, TIM50
```

### Coordinate conversion

The script converts smallFISH `[z, y, x]` pixel coordinates to `[x_um, y_um, z_um]` using:

```python
X_SCALE = 0.0645
Y_SCALE = 0.0645
Z_SCALE = 0.2
```

### Output per `spots_extraction` folder

```text
converted_coordinates/
    cell_000_xyz_um.csv
    cell_001_xyz_um.csv
    ...

<mRNA>_NN_distance_um.npy
<mRNA>_NN_distance_um.xlsx
<mRNA>_spot_counts.csv
<mRNA>_spot_counts.xlsx
<mRNA>_NN_distribution.png
```

If no nearest-neighbor distances are possible, for example if all cells have fewer than two spots, it writes:

```text
<mRNA>_NO_NN_DISTANCES.txt
```

### Parent-level output

```text
Condition parent folder/NN_distance_recursive_run_summary.csv
Condition parent folder/NN_distance_recursive_run_summary.xlsx
```

### Important downstream files

- `converted_coordinates/cell_###_xyz_um.csv`  
  Real mRNA coordinates in microns. Used by optional node-distance scripts.

- `<mRNA>_NN_distance_um.npy` and `<mRNA>_NN_distance_um.xlsx`  
  Real mRNA nearest-neighbor distance distribution. The `.npy` file is used by Python downstream scripts; the `.xlsx` file is for easy inspection in Excel.

- `<mRNA>_spot_counts.csv` and `<mRNA>_spot_counts.xlsx`  
  Per-cell mRNA counts. The CSV is used by random null scripts to generate equal-count random points; the Excel version is for easy inspection.

- `<mRNA>_NN_distribution.png`  
  QC histogram of real mRNA nearest-neighbor distances.

---


## 12A. `rna_intensity_histograms_recursive.py`

### Purpose

Creates RNA spot-intensity histograms from raw smallFISH `spots_extractions` files. It summarizes intensity distributions by Series, by experimental condition, and by RNA channel.

This script is useful for detecting intensity differences caused by acquisition settings, thresholding, probe quality, or sample quality. It is especially useful when comparing RNA detection across multiple Series or experimental conditions.

### Input

```text
Condition or experiment parent folder/
└── Series XX/
    └── <mRNA>/.../results/spots_extraction/
        └── spots_extractions_*.csv or .xlsx
```

The script searches recursively for files whose names contain:

```text
spots_extractions
```

It auto-detects the intensity column. The usual smallFISH column is:

```text
intensity
```

The script skips generated downstream folders such as:

```text
converted_coordinates
random_*_output
plots
pooled_real_vs_random_comparisons
mrna_coloc_BY_SERIES
mito_rna_surface_GLOBAL_CALIBRATION
node_mrna_dist_BY_SERIES
```

This helps prevent generated output tables from being accidentally reprocessed as raw smallFISH input.

### Key settings

```python
HISTOGRAM_BINS = 50
INTENSITY_XMIN = None
INTENSITY_XMAX = None
MAKE_LOG10_HISTOGRAMS = True
WRITE_EXCEL = True
MAKE_SERIES_POOLED_HISTOGRAMS = True
MAKE_SERIES_BY_CHANNEL_HISTOGRAMS = True
MAKE_CONDITION_POOLED_HISTOGRAMS = True
MAKE_CONDITION_BY_CHANNEL_HISTOGRAMS = True
MAKE_ALL_PARENT_HISTOGRAM = True
```

Optional intensity filtering settings:

```python
FILTER_MIN_INTENSITY = None
FILTER_MAX_INTENSITY = None
```

Keep these as `None` unless there is a documented reason to exclude intensity outliers.

### Output folder

```text
RNA_intensity_histograms/
```

### Main workbook

```text
RNA_intensity_histogram_analysis_workbook.xlsx
```

This workbook combines the major summary tables:

| Sheet | Meaning |
|---|---|
| `all_intensities` | One row per RNA spot, with condition, Series, channel, source file, cell index, and intensity. |
| `file_summary` | One row per raw smallFISH file, including number of finite intensities extracted. |
| `series_summary` | Summary statistics by condition and Series. |
| `series_channel` | Summary statistics by condition, Series, and RNA channel. |
| `condition_summary` | Summary statistics by condition. |
| `condition_channel` | Summary statistics by condition and RNA channel. |
| `channel_summary` | Summary statistics by RNA channel across the selected parent folder. |
| `rejected_files` | Files that were rejected or failed during parsing. |
| `settings` | Exact run settings. |

### Individual table outputs

The same information is also saved as separate CSV and Excel files:

```text
ALL_RNA_spot_intensities.csv/.xlsx
RNA_intensity_file_summary.csv/.xlsx
RNA_intensity_summary_by_series.csv/.xlsx
RNA_intensity_summary_by_series_and_channel.csv/.xlsx
RNA_intensity_summary_by_condition.csv/.xlsx
RNA_intensity_summary_by_condition_and_channel.csv/.xlsx
RNA_intensity_summary_by_channel.csv/.xlsx
RNA_intensity_plot_index.csv/.xlsx
RNA_intensity_rejected_files.csv/.xlsx
RNA_intensity_RUN_SETTINGS.csv/.xlsx
```

### Plot outputs

```text
plots/parent/
plots/by_series/
plots/by_series_and_channel/
plots/by_condition/
plots/by_condition_and_channel/
```

Important PNGs include:

```text
ALL_PARENT_RNA_intensity_histogram.png
ALL_PARENT_RNA_intensity_histogram_log10.png
RNA_intensity_by_series__condition_label_<condition>__series_label_<series>.png
RNA_intensity_by_condition__condition_label_<condition>.png
RNA_intensity_by_condition_and_channel__condition_label_<condition>__rna_channel_<channel>.png
```

### How to interpret

Use the summary tables to compare mean, median, standard deviation, SEM, quartiles, and min/max intensity values across Series and conditions. Use the histograms to identify shifts, long tails, bimodal distributions, or outlier-heavy Series.

If one condition has systematically higher or lower intensity distributions, check whether imaging settings, thresholding, exposure, sample preparation, or probe performance differed before interpreting spatial differences as biological effects.

---


## 13. `mrna_colocalization_by_cell.py`

### Purpose

Compares two mRNA channels cell-by-cell within each Series folder. It answers questions such as:

```text
Which MS2 puncta overlap with ATP3 puncta within 0.25 um?
What percent of MS2 puncta have a nearby ATP3 punctum in the same cell?
What percent of ATP3 puncta have a nearby MS2 punctum in the same cell?
```

The script uses a 3D Euclidean distance cutoff. A punctum is called colocalized if its nearest punctum from the other channel is within the adjustable threshold.

### Key settings

```python
MRNA_A = "MS2"
MRNA_B = "ATP3"
COLOCALIZATION_DISTANCE_UM = 0.25
```

Change `MRNA_A`, `MRNA_B`, and `COLOCALIZATION_DISTANCE_UM` to compare a different pair or distance range. Common test thresholds are `0.10`, `0.25`, `0.50`, and `1.00` microns.

### Input

This script uses converted coordinate files from `6.mrna_nn_distance_recursive.py`:

```text
Condition parent folder/
└── Series XX/
    ├── MS2/.../results/spots_extraction/converted_coordinates/cell_###_xyz_um.csv
    └── ATP3/.../results/spots_extraction/converted_coordinates/cell_###_xyz_um.csv
```

The matching is strict by Series and cell index:

```text
MS2 cell_000 compares only to ATP3 cell_000 in the same Series folder.
MS2 cell_001 compares only to ATP3 cell_001 in the same Series folder.
```

### Output folder

For the default settings, the main output folder is:

```text
Condition parent folder/mrna_coloc_BY_SERIES/MS2_vs_ATP3_within_0p25um/
```

The distance tag changes with the selected cutoff. For example, a `0.50 um` cutoff produces:

```text
MS2_vs_ATP3_within_0p5um/
```

### Main all-cell outputs

```text
ALL_cell_colocalization_summary.csv
ALL_cell_colocalization_summary.xlsx
ALL_MS2_nearest_ATP3.csv
ALL_MS2_nearest_ATP3.xlsx
ALL_ATP3_nearest_MS2.csv
ALL_ATP3_nearest_MS2.xlsx
ALL_MS2_ATP3_pairs_within_0p25um.csv
ALL_MS2_ATP3_pairs_within_0p25um.xlsx
ALL_MS2_vs_ATP3_colocalization_workbook.xlsx
ALL_unmatched_or_failed_cells.csv
ALL_unmatched_or_failed_cells.xlsx
ALL_duplicate_coordinate_files.csv
ALL_duplicate_coordinate_files.xlsx
ALL_input_coordinate_file_index.csv
ALL_input_coordinate_file_index.xlsx
RUN_SETTINGS.csv
RUN_SETTINGS.xlsx
```

### Most important workbook

```text
ALL_MS2_vs_ATP3_colocalization_workbook.xlsx
```

This workbook contains multiple sheets:

| Sheet | Meaning |
|---|---|
| `cell_summary` | One row per matched cell with counts, percent colocalized, and summary distances. |
| `MS2_nearest_ATP3` | One row per MS2 punctum; gives the nearest ATP3 punctum and whether it is within the cutoff. |
| `ATP3_nearest_MS2` | One row per ATP3 punctum; gives the nearest MS2 punctum and whether it is within the cutoff. |
| `pairs_within_threshold` | Every MS2-ATP3 punctum pair whose 3D distance is within the selected cutoff. |
| `unmatched_failed` | Cells that could not be analyzed because one channel was missing or a file failed. |
| `duplicates` | Duplicate coordinate files for the same Series, mRNA, and cell index. |
| `input_index` | Audit table of all coordinate files found by the script. |
| `settings` | Exact run settings, including mRNA names, cutoff, scales, parent folder, and output folder. |

### `ALL_cell_colocalization_summary.xlsx`

This is the best file for downstream statistics. It has one row per matched cell.

Important columns:

| Column | Meaning |
|---|---|
| `series_label` | Series folder analyzed. |
| `cell_index` / `cell_label` | Cell number, such as `0` / `cell_000`. |
| `colocalization_distance_um` | Distance cutoff used for the run. |
| `A_count` | Number of puncta in `MRNA_A`, for example MS2. |
| `B_count` | Number of puncta in `MRNA_B`, for example ATP3. |
| `A_colocalized_count_nearest_neighbor` | Number of A puncta whose nearest B punctum is within the cutoff. |
| `B_colocalized_count_nearest_neighbor` | Number of B puncta whose nearest A punctum is within the cutoff. |
| `A_percent_colocalized_nearest_neighbor` | Percent of A puncta colocalized with B. |
| `B_percent_colocalized_nearest_neighbor` | Percent of B puncta colocalized with A. |
| `pair_count_within_threshold` | Total number of A-B punctum pairs within the cutoff. |
| `unique_A_puncta_in_any_pair` | Number of unique A puncta involved in at least one close pair. |
| `unique_B_puncta_in_any_pair` | Number of unique B puncta involved in at least one close pair. |
| `A_to_B_median_nearest_distance_um` | Median nearest-B distance for A puncta. |
| `B_to_A_median_nearest_distance_um` | Median nearest-A distance for B puncta. |

Nearest-neighbor colocalization is directional. The percent of MS2 near ATP3 does not have to equal the percent of ATP3 near MS2, because one ATP3 punctum can be the nearest neighbor for multiple MS2 puncta.

### `ALL_MS2_nearest_ATP3.xlsx`

This file has one row per MS2 punctum. It reports the nearest ATP3 punctum in the same cell.

Important columns:

```text
MS2_index
MS2_x, MS2_y, MS2_z
nearest_ATP3_index
nearest_ATP3_x, nearest_ATP3_y, nearest_ATP3_z
nearest_ATP3_distance_um
colocalized
```

Use this file to ask: for every MS2 punctum, how close is the nearest ATP3 punctum?

### `ALL_ATP3_nearest_MS2.xlsx`

This is the reverse nearest-neighbor table. It has one row per ATP3 punctum and reports the nearest MS2 punctum in the same cell.

Use this file to ask: for every ATP3 punctum, how close is the nearest MS2 punctum?

### `ALL_MS2_ATP3_pairs_within_0p25um.xlsx`

This file lists every exact MS2-ATP3 punctum pair whose distance is within the cutoff. It can contain more rows than the number of MS2 or ATP3 puncta because one punctum can pair with multiple puncta from the other channel.

Important columns:

```text
MS2_index
MS2_x, MS2_y, MS2_z
ATP3_index
ATP3_x, ATP3_y, ATP3_z
pair_distance_um
```

Use this file when you need the exact identities and distances of overlapping puncta.

### QC and audit outputs

| File | Meaning |
|---|---|
| `ALL_unmatched_or_failed_cells.xlsx` | Cells where one channel was missing or a file failed to read. Ideally, this is empty or explainable. |
| `ALL_duplicate_coordinate_files.xlsx` | Duplicate coordinate files for the same Series/mRNA/cell. Check this if outputs look duplicated. |
| `ALL_input_coordinate_file_index.xlsx` | Every coordinate file found and how it was assigned to Series, mRNA, and cell index. |
| `RUN_SETTINGS.xlsx` | Exact settings used for reproducibility. |

### Plots

```text
ALL_percent_colocalized_by_cell.png
ALL_MS2_nearest_ATP3_distance_histogram.png
ALL_ATP3_nearest_MS2_distance_histogram.png
ALL_MS2_ATP3_pair_distances_within_0p25um.png
```

The per-cell plot folder is:

```text
per_series/<SeriesLabel>/cell_###/
```

Each cell folder contains:

```text
cell_###_MS2_vs_ATP3_MS2_nearest_ATP3.xlsx
cell_###_MS2_vs_ATP3_ATP3_nearest_MS2.xlsx
cell_###_MS2_vs_ATP3_pairs_within_0p25um.xlsx
cell_###_MS2_vs_ATP3_summary.xlsx
cell_###_MS2_vs_ATP3_xy_colocalization.png
```

The XY plot is only a 2D projection for visual QC. The colocalization call itself uses 3D x, y, z distance.

## 14. `7.random_mrna_null_recursive.py`

### Purpose

Generates the standard random mRNA null model. It keeps the real mRNA count per cell but randomizes mRNA positions near the MitoGraph skeleton from the same series.

### Randomization rule

For each real mRNA count table row:

1. Match the cell to a MitoGraph skeleton `.txt` file in the same `Series XX/cells` folder.
2. Load skeleton `x, y, z` coordinates.
3. For each mRNA count, choose a random skeleton coordinate.
4. Choose a random 3D direction.
5. Place the random point within:

```python
MAX_DISTANCE_UM = 0.5
```

of the selected skeleton coordinate.

6. Calculate nearest-neighbor distances among the generated random mRNA points.

### Input

```text
Series XX/cells/*.txt
Series XX/<mRNA>/.../results/spots_extraction/<mRNA>_spot_counts.csv
```

The script builds one shared skeleton index per `Series XX` folder so both MS2 and ATP2/ATP3/TIM50 use the same mitochondrial skeleton files for that series.

### Output per `spots_extraction` folder

```text
random_<mRNA>_output/
    cell_000_random_<mRNA>.csv
    cell_001_random_<mRNA>.csv
    ...
    random_<mRNA>_nn.npy
    random_<mRNA>_nn.csv
    random_<mRNA>_distribution.png
    random_<mRNA>_per_cell_summary.csv
```

If no nearest-neighbor distances are generated:

```text
random_<mRNA>_NO_NN_DISTANCES.txt
```

### Parent-level output

```text
Condition parent folder/random_mrna_recursive_run_summary_SHARED_SKELETON.csv
```

### Important downstream file

```text
random_<mRNA>_output/random_<mRNA>_nn.npy
```

This is the random distribution expected by `8.real_vs_random_recursive_comparison.py`, `9.pooled_real_vs_random_across_series_REPLICATE_LEVEL_STATS.py`, and `10. pooled_real_vs_random_ Thesis Graphs.py`.

---

## 15. `7.1 random_mrna_distance_range.py`

### Purpose

Alternative random null-model script for sensitivity analyses. It allows a configurable distance range from skeleton coordinates instead of only `0` to `MAX_DISTANCE_UM`.

### Key settings

```python
MIN_DISTANCE_UM = 0.0
MAX_DISTANCE_UM = 5.0
ENFORCE_NEAREST_SKELETON_RANGE = False
```

Examples:

```python
MIN_DISTANCE_UM = 0.0
MAX_DISTANCE_UM = 5.0
```

places random points 0 to 5 microns from the chosen skeleton coordinate.

```python
MIN_DISTANCE_UM = 5.0
MAX_DISTANCE_UM = 5.0
```

places random points exactly 5 microns from the chosen skeleton coordinate.

If `ENFORCE_NEAREST_SKELETON_RANGE = True`, the script also checks the nearest skeleton coordinate after generation and rejects candidates outside the requested distance range. This is stricter but slower.

### Input

Same as `7.random_mrna_null_recursive.py`:

```text
Series XX/cells/*.txt
Series XX/<mRNA>/.../results/spots_extraction/<mRNA>_spot_counts.csv
```

### Output per `spots_extraction` folder

The output folder includes a distance tag:

```text
random_<mRNA>_output_dist_<MIN>_to_<MAX>um/
```

Example:

```text
random_ATP3_output_dist_0_to_5um/
random_MS2_output_dist_5_to_5um/
```

Files inside the folder also include the distance tag:

```text
cell_000_random_<mRNA>_dist_<MIN>_to_<MAX>um.csv
random_<mRNA>_nn_dist_<MIN>_to_<MAX>um.npy
random_<mRNA>_nn_dist_<MIN>_to_<MAX>um.csv
random_<mRNA>_distribution_dist_<MIN>_to_<MAX>um.png
random_<mRNA>_per_cell_summary_dist_<MIN>_to_<MAX>um.csv
```

### Parent-level output

```text
Condition parent folder/random_mrna_recursive_run_summary_SHARED_SKELETON_dist_<MIN>_to_<MAX>um.csv
```

### Important compatibility warning

The standard downstream comparison scripts look for:

```text
random_<mRNA>_output/random_<mRNA>_nn.npy
```

The distance-range script writes distance-tagged files instead. Use it for sensitivity analysis unless the downstream comparison script has been modified to find the distance-tagged null file. If you want a distance-range output to be treated as the default random distribution, copy or rename the selected distance-tagged `.npy` into the standard expected location only after documenting which distance range was used.

---

## 16. `8.real_vs_random_recursive_comparison.py`

### Purpose

Makes a real-vs-random comparison for each mRNA and each `spots_extraction` folder.

This is a per-series/per-folder QC step. It helps identify whether a particular series has abnormal behavior before pooled condition-level analysis.

### Input

For each mRNA-specific `spots_extraction` folder:

```text
<mRNA>_NN_distance_um.npy
random_<mRNA>_output/random_<mRNA>_nn.npy
```

Examples:

```text
MS2_NN_distance_um.npy
random_MS2_output/random_MS2_nn.npy

ATP3_NN_distance_um.npy
random_ATP3_output/random_ATP3_nn.npy
```

### Output per `spots_extraction` folder

```text
plots/
    <mRNA>_real_vs_random.png
    <mRNA>_real_vs_random_summary.csv
```

If the random file is missing:

```text
plots/<mRNA>_NO_RANDOM_FILE_FOUND.txt
```

If the real or random array is empty:

```text
plots/<mRNA>_real_vs_random_NOT_CREATED.txt
```

### Parent-level output

```text
Condition parent folder/real_vs_random_recursive_comparison_summary.csv
```

### Output interpretation

The PNG overlays density-normalized real and random nearest-neighbor distance histograms for that mRNA. The summary CSV stores basic counts, means, medians, and standard deviations for that single comparison.

---

## 17. `9.pooled_real_vs_random_across_series_REPLICATE_LEVEL_STATS.py`

### Purpose

Pools all valid series within one condition folder and creates one pooled real-vs-random comparison per mRNA.

This is the preferred condition-level statistical analysis because it uses replicate/series-level tests instead of treating every mRNA spot as an independent biological replicate.

### Input

The selected parent folder should contain all `Series XX` folders for one condition. The script searches recursively for real files:

```text
*_NN_distance_um.npy
```

and matches them to standard random files:

```text
random_<mRNA>_output/random_<mRNA>_nn.npy
```

### Key statistical settings

```python
DIFFRACTION_LIMIT_UM = 0.20
PRACTICAL_EFFECT_THRESHOLD_UM = 0.05
MAIN_PLOT_TEST = "replicate_median_permutation"
N_PERMUTATIONS = 20000
```

Distances below `0.20 µm` are floored to `0.20 µm` for plotted and statistical diffraction-limited values.

### Output folder

```text
Condition parent folder/pooled_real_vs_random_comparisons/
```

### Output files per mRNA

```text
<mRNA>_POOLED_real_NN_distance_um_RAW.npy
<mRNA>_POOLED_random_NN_distance_um_RAW.npy
<mRNA>_POOLED_real_NN_distance_um_DIFFRACTION_FLOORED.npy
<mRNA>_POOLED_random_NN_distance_um_DIFFRACTION_FLOORED.npy
<mRNA>_POOLED_real_vs_random_values.csv
<mRNA>_REPLICATE_LEVEL_real_vs_random_summary.csv
<mRNA>_POOLED_real_vs_random_REPLICATE_STATS.png
```

### Output files for the condition

```text
pooled_input_file_pairing_summary.csv
pooled_real_vs_random_summary.csv
```

### How to interpret

Use the replicate-level permutation p-value and replicate-level effect columns for biological interpretation. The point-level Mann-Whitney and KS statistics are saved in the CSV, but they can be overly sensitive when many mRNA spots are pooled.

---

## 18. `10. pooled_real_vs_random_ Thesis Graphs.py`

### Purpose

Creates thesis-ready pooled real-vs-random plots. It is based on the pooled replicate-level script but adds cleaner figure styling, strain-labeled CSV filenames, median-line annotations, and a stricter inference gate.

### Input

Same as script 9:

```text
Condition parent folder containing Series folders
*_NN_distance_um.npy
random_<mRNA>_output/random_<mRNA>_nn.npy
```

### Key settings

```python
DIFFRACTION_LIMIT_UM = 0.20
PRACTICAL_EFFECT_THRESHOLD_UM = 0.10
STRICT_ALPHA = 0.001
STRICT_MIN_REPLICATE_PAIRS = 3
STRICT_MIN_DIRECTION_AGREEMENT = 0.75
STRICT_REQUIRE_CI_EXCLUDES_ZERO = True
STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD = True
MAIN_PLOT_TEST = "replicate_median_permutation"
N_PERMUTATIONS = 20000
```

### Output folder

```text
Condition parent folder/pooled_real_vs_random_comparisons/
```

### Output files per mRNA

```text
<mRNA>_POOLED_real_NN_distance_um_RAW.npy
<mRNA>_POOLED_random_NN_distance_um_RAW.npy
<mRNA>_POOLED_real_NN_distance_um_DIFFRACTION_FLOORED.npy
<mRNA>_POOLED_random_NN_distance_um_DIFFRACTION_FLOORED.npy
<mRNA>_POOLED_real_vs_random_values_<strain>.csv
<mRNA>_REPLICATE_LEVEL_real_vs_random_summary_<strain>.csv
<mRNA>_POOLED_real_vs_random_REPLICATE_STATS_<strain>.png
```

### Output files for the condition

```text
pooled_input_file_pairing_summary_<strain>.csv
pooled_real_vs_random_summary_<strain>.csv
```

### How to interpret

Use this version for final figure generation. The plot itself is intentionally cleaner. Detailed p-values, confidence intervals, strict-call status, strict-call reason, and replicate-level effect metrics are saved in the CSV files.

---

## 19. `11.compare_strains_probe_sets_boxplot_v6.py`

### Purpose

Compares strain/probe effects across many condition folders. This is the cross-condition summary step.

The script expects each strain/probe folder to already contain a `pooled_real_vs_random_comparisons` folder with a pairing summary. It combines series-level real-vs-random metrics into grouped boxplots and statistical summary tables.

### Expected parent layout

Select the folder that contains strain folders:

```text
Experiment parent folder/
├── yWL333/
│   ├── MS2(ATP6)/pooled_real_vs_random_comparisons/
│   ├── MS2(ATP8)/pooled_real_vs_random_comparisons/
│   ├── MS2(ATP2)/pooled_real_vs_random_comparisons/
│   ├── MS2(ATP3)/pooled_real_vs_random_comparisons/
│   └── MS2(TIM50)/pooled_real_vs_random_comparisons/
├── yMM002(ATP11)/
│   ├── MS2(ATP6_8)/pooled_real_vs_random_comparisons/
│   ├── ATP2/pooled_real_vs_random_comparisons/
│   ├── ATP3/pooled_real_vs_random_comparisons/
│   └── TIM50/pooled_real_vs_random_comparisons/
└── ...
```

The script also supports repeated inner strain folders such as:

```text
Experiment parent/strain/probe/strain/pooled_real_vs_random_comparisons/
```

### Required input inside each pooled folder

```text
pooled_input_file_pairing_summary.csv
```

or equivalent strain-suffixed pairing files if the script is pointed at folders created by the thesis-graph version and the file naming is adjusted accordingly.

### Probe assignment behavior

ATP6 and ATP8 are collapsed into one probe set:

```text
ATP6/8
```

Rows labeled `MS2` are interpreted as the ATP6/8 reporter. In mixed folders such as `MS2 (ATP2)`, the script first checks the row-level `mRNA` label, so MS2 rows are assigned to ATP6/8 while ATP2 rows remain ATP2.

### Output folder

```text
Experiment parent folder/strain_probe_comparison/
```

### Main output files

```text
strain_probe_folder_mapping.csv
strain_probe_series_level_metrics.csv
probe_set_mRNA_assignment_counts.csv
strain_probe_summary.csv
statistics_by_probe_across_strains_delta_median.csv
statistics_by_strain_across_probes_delta_median.csv
wt_vs_ko_statistics_delta_median.csv
strain_probe_delta_median_grouped_boxplot.png
```

If enabled, SVG copies are also saved:

```text
strain_probe_delta_median_grouped_boxplot.svg
```

Individual probe plots are saved when `make_individual_probe_plots = True`:

```text
ATP6_8_delta_median_by_strain_boxplot.png
ATP2_delta_median_by_strain_boxplot.png
ATP3_delta_median_by_strain_boxplot.png
TIM50_delta_median_by_strain_boxplot.png
```

Optional outputs, depending on settings:

```text
strain_probe_delta_median_grouped_by_strain_boxplot.png
strain_probe_delta_median_summary_heatmap.png
```

### Main metric

Default:

```python
primary_metric = "delta_median"
```

where:

```text
delta_median = median(real) - median(random)
```

The grouped boxplot uses one dot per series and box/whisker summaries per strain within each probe set.

---

# PART E — Optional mitochondrial proximity and node-distance scripts

These scripts answer a different question from the mRNA nearest-neighbor/random-null pipeline. They calculate how far each mRNA spot is from the nearest MitoGraph node or node component.

Use these scripts when you want mRNA-to-mitochondrial-node proximity, not mRNA-to-mRNA nearest-neighbor spacing.

---


## 20. `mito_rna_surface_GLOBAL_CALIBRATION.py`

### Purpose

Measures how close each mRNA spot is to the mitochondrial surface using one fixed global calibration. This answers a different question from mRNA-mRNA colocalization and mRNA nearest-neighbor spacing:

```text
How close are MS2, ATP2, ATP3, or Tim50 puncta to the mitochondrial surface?
What fraction of spots are within 0.25, 0.5, 0.75, or 1.0 um of the mito surface?
```

The script also creates interactive 3D HTML overlays of the mitochondrial surface, skeleton, RNA spots, and nearest-surface links.

### Important calibration behavior

This script does not perform per-cell fitting. It applies the same fixed transform to all cells and channels:

```python
FIXED_TRANSFORM = "flip_y"
X_SCALE = 0.05805
Y_SCALE = 0.05805
Z_SCALE = 0.2
GLOBAL_DX_UM = -0.125
GLOBAL_DY_UM = -0.250
GLOBAL_DZ_UM = 0.000
CROP_WIDTH_PIXELS = 200
CROP_HEIGHT_PIXELS = 200
```

The fixed calibration came from the MS2-positive-control scale/offset diagnostic and is applied identically to MS2, ATP2, ATP3, and Tim50.

### Channels and thresholds

Default channels:

```python
RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "Tim50"]
```

Default mitochondrial-surface thresholds:

```python
COLOCALIZATION_THRESHOLDS_UM = [0.25, 0.5, 0.75, 1.0]
```

These thresholds are used to summarize the fraction and number of RNA spots close to the mitochondrial surface.

### Input

This script reads raw smallFISH `spots_extractions` files and MitoGraph mitochondrial surface files.

Expected inputs:

```text
Condition parent folder/
└── Series XX/
    ├── cells/
    │   ├── *_mitosurface.vtk
    │   └── MitoGraph skeleton .txt files
    ├── MS2/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    ├── ATP2/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    ├── ATP3/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    └── Tim50/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
```

The script intentionally skips generated folders such as:

```text
converted_coordinates
random_*_output
mito_visualization_interactive
mito_rna_surface_GLOBAL_CALIBRATION
```

This prevents generated outputs from being reprocessed as raw smallFISH inputs.

### Output folder

At the condition parent level:

```text
Condition parent folder/mito_rna_surface_GLOBAL_CALIBRATION/
```

Inside each `cells` folder:

```text
Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/
```

### Parent-level outputs

```text
PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv
PARENT_spot_level_surface_distances_GLOBAL_CALIBRATION.csv
PARENT_channel_summary_GLOBAL_CALIBRATION.csv
PARENT_cell_level_surface_colocalization_CORRECTED_XY_0p065_CROP250.csv
PARENT_spot_level_surface_distances_CORRECTED_XY_0p065_CROP250.csv
PARENT_channel_summary_CORRECTED_XY_0p065_CROP250.csv
CORRECTED_XY_0p065_CROP250_RUN_SETTINGS.csv
```

### `PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv`

This is one of the main summary files. It has one row per Series/cell/channel combination.

Important columns include:

| Column | Meaning |
|---|---|
| `SeriesFolder` | Series folder analyzed. |
| `CellsFolder` | `cells` folder containing MitoGraph surface files. |
| `CellIndex` | Matched cell index. |
| `Channel` | RNA channel, such as MS2, ATP2, ATP3, or Tim50. |
| `MitoSurfaceFile` | MitoGraph `_mitosurface.vtk` used for that cell. |
| `RNAFile` | Raw smallFISH spot file used for that cell/channel. |
| `Transform` | Fixed transform, usually `flip_y`. |
| `X_SCALE_um_per_px`, `Y_SCALE_um_per_px`, `Z_SCALE_um_per_px` | Pixel-to-micron scale values. |
| `GlobalDX_um`, `GlobalDY_um`, `GlobalDZ_um` | Global offset applied after scaling/flipping. |
| `Status` | Processing status, usually `processed_global_calibration`. |
| `Nspots` | Number of RNA spots in that cell/channel. |
| `MeanDistanceToSurface_um` | Mean RNA-to-surface distance. |
| `MedianDistanceToSurface_um` | Median RNA-to-surface distance. |
| `MinDistanceToSurface_um` | Closest RNA-to-surface distance. |
| `MaxDistanceToSurface_um` | Largest RNA-to-surface distance. |
| `FractionWithin0p25um`, `SpotsWithin0p25um` | Fraction/count within 0.25 um of the surface. |
| `FractionWithin0p5um`, `SpotsWithin0p5um` | Fraction/count within 0.5 um of the surface. |
| `FractionWithin0p75um`, `SpotsWithin0p75um` | Fraction/count within 0.75 um of the surface. |
| `FractionWithin1p0um`, `SpotsWithin1p0um` | Fraction/count within 1.0 um of the surface. |

### `PARENT_spot_level_surface_distances_GLOBAL_CALIBRATION.csv`

This is the spot-level file. It contains one row per RNA punctum.

Important columns include:

```text
x_um, y_um, z_um
x_px, y_px, z_px
DistanceToMitoSurface_um
nearest_surface_x_um, nearest_surface_y_um, nearest_surface_z_um
within_0p25_um
within_0p5_um
within_0p75_um
within_1p0_um
SeriesFolder
CellsFolder
CellIndex
Channel
MitoSurfaceFile
RNAFile
```

Use this file when you need individual RNA spot distances to the mitochondrial surface.

### `PARENT_channel_summary_GLOBAL_CALIBRATION.csv`

This aggregates across cells by RNA channel. It reports cell counts, total spots, and mean/median/SEM across cells for the distance and threshold metrics.

Use this file when comparing overall mitochondrial-surface proximity across channels.

### Per-cell outputs inside each `cells` folder

```text
GLOBAL_CALIBRATION_ACCEPTED_raw_spot_files.csv
GLOBAL_CALIBRATION_REJECTED_nonraw_files.csv
cell_###_<Channel>_GLOBAL_CALIBRATION_surface_distances.csv
cell_###_GLOBAL_CALIBRATION_surface_overlay.html
```

Output meanings:

| File | Meaning |
|---|---|
| `GLOBAL_CALIBRATION_ACCEPTED_raw_spot_files.csv` | Raw smallFISH files accepted for analysis. |
| `GLOBAL_CALIBRATION_REJECTED_nonraw_files.csv` | Candidate files rejected because they were generated outputs or not raw spot files. |
| `cell_###_<Channel>_GLOBAL_CALIBRATION_surface_distances.csv` | Per-spot transformed coordinates and distances to the mito surface for one cell/channel. |
| `cell_###_GLOBAL_CALIBRATION_surface_overlay.html` | Interactive 3D overlay showing mitochondrial surface, optional skeleton, RNA spots, and nearest-surface links. |

### Output interpretation

The distance is computed to the nearest mitochondrial surface vertex from the MitoGraph `_mitosurface.vtk` mesh. This is an approximation of distance to the surface mesh and is best interpreted as a consistent proximity metric across cells/channels, not a manual colocalization score.


## 20A. `mito_rna_surface_CORRECTED_XY_0p065_CROP250.py`

### Purpose

This is the current corrected version of the mitochondrial-surface proximity script for datasets where MitoGraph was processed with the corrected microscope calibration. Use this version when the old `GLOBAL_CALIBRATION` surface overlays are visibly shifted or when the MitoGraph command used the corrected XY calibration.

It measures the 3D distance from each raw smallFISH RNA spot to the nearest mitochondrial surface point and creates interactive 3D overlays for visual QC.

### When to use this version

Use this version instead of the older `mito_rna_surface_GLOBAL_CALIBRATION.py` when:

- MitoGraph was run with corrected microscope XY calibration,
- surface overlays from the older script look shifted,
- RNA points appear too low after `flip_y`,
- the cropped cell images are 250 x 250 pixels,
- you want corrected output files kept separate from older global-calibration outputs.

### Corrected calibration settings

Default settings:

```python
CORRECTED_XY_SCALE_UM_PER_PIXEL = 0.065
CORRECTED_Z_SCALE_UM_PER_PIXEL = 0.2
FIXED_TRANSFORM = "flip_y"
X_SCALE = 0.065
Y_SCALE = 0.065
Z_SCALE = 0.2
GLOBAL_DX_UM = 0.000
GLOBAL_DY_UM = 0.000
GLOBAL_DZ_UM = 0.000
CROP_WIDTH_PIXELS = 250
CROP_HEIGHT_PIXELS = 250
```

The crop-size correction is important because `flip_y` uses the crop height:

```python
y_um = CROP_HEIGHT_PIXELS * Y_SCALE - y_px * Y_SCALE
```

Changing the crop height from 200 to 250 pixels shifts flipped-Y RNA coordinates upward by:

```text
(250 - 200) x 0.065 = 3.25 um
```

That large shift is why this corrected version should be used when the RNA overlay appears far below the mitochondrial surface.

### Transform options

The default transform is:

```python
FIXED_TRANSFORM = "flip_y"
```

If the scale is correct but the overlay is mirrored, test only one transform change at a time:

```python
FIXED_TRANSFORM = "none"
FIXED_TRANSFORM = "flip_x"
FIXED_TRANSFORM = "flip_xy"
```

Do not change scale, crop size, and transform simultaneously. Change one setting, rerun one representative cell or Series, and inspect the HTML overlay.

### Input

Expected inputs are the same as the older surface script:

```text
Condition parent folder/
└── Series XX/
    ├── cells/
    │   ├── *_mitosurface.vtk
    │   └── MitoGraph skeleton .txt files
    ├── MS2/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    ├── ATP2/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    ├── ATP3/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
    └── Tim50/.../results/spots_extraction/spots_extractions_*.csv or .xlsx
```

### Output folder

At the condition parent level:

```text
Condition parent folder/mito_rna_surface_CORRECTED_XY_0p065_CROP250/
```

Inside each `cells` folder:

```text
Series XX/cells/mito_rna_surface_CORRECTED_XY_0p065_CROP250/
```

### Main outputs

```text
PARENT_cell_level_surface_colocalization_CORRECTED_XY_0p065_CROP250.csv
PARENT_spot_level_surface_distances_CORRECTED_XY_0p065_CROP250.csv
PARENT_channel_summary_CORRECTED_XY_0p065_CROP250.csv
CORRECTED_XY_0p065_CROP250_RUN_SETTINGS.csv
```

Per-cell outputs include:

```text
cell_###_<Channel>_CORRECTED_XY_0p065_CROP250_surface_distances.csv
cell_###_CORRECTED_XY_0p065_CROP250_surface_overlay.html
```

### Output interpretation

Use the CSV files for quantitative analysis and the HTML overlays for alignment QC. The overlays are especially important for this script because a wrong crop height, mirror transform, or offset can produce plausible-looking distances that are biologically wrong.

Check at least a few cells per condition before trusting the surface-distance numbers.

---

## 21. `mrna_to_nearest_node_distance_BY_SERIES_FIXED.py`

### Purpose

Computes distances from real mRNA coordinates to the nearest MitoGraph node for each series. This fixed version is the recommended node-distance script.

### Why the fixed version matters

The earlier node-distance script could misread the `-002-1` part of an ND2 filename as a cell index and accidentally match `cell_002` to `_000_nodes.vtk`. The fixed script extracts the true cell index from the final `_000_nodes.vtk` suffix and matches:

```text
cell_000_xyz_um.csv -> *_000_nodes.vtk
cell_001_xyz_um.csv -> *_001_nodes.vtk
cell_002_xyz_um.csv -> *_002_nodes.vtk
```

It also avoids Windows path-length errors by using shorter output folder names and excludes derived surface-distance tables.

### Input

```text
Condition parent folder/
└── Series XX/
    ├── cells/*_nodes.vtk
    └── <mRNA>/.../results/spots_extraction/converted_coordinates/cell_*_xyz_um.csv
```

By default, it only processes raw converted-coordinate files:

```text
converted_coordinates/cell_*_xyz_um.csv
```

It excludes derived tables such as:

```text
mito_rna_surface_GLOBAL_CALIBRATION/*_surface_distances.csv
```

### Key setting

```python
NODE_POSITION_MODE = "connected_component_centers"
```

This treats each connected node mesh component as one node coordinate by using its center.

### Output folder

```text
Condition parent folder/node_mrna_dist_BY_SERIES/
```

### Per-series outputs

```text
node_mrna_dist_BY_SERIES/per_series/<SeriesLabel>/
    <SeriesLabel>_node_file_index.csv
    <SeriesLabel>_coordinate_file_index.csv
    <mRNA>_cell_###/
        <mRNA>_cell_###_distances.csv
        <mRNA>_cell_###_node_occupancy.csv
        <mRNA>_cell_###_summary.csv
        <mRNA>_cell_###_histogram.png
        <mRNA>_cell_###_threshold_counts.png
        <mRNA>_cell_###_xy_distance.png
        <mRNA>_cell_###_node_occupancy.png
    <SeriesLabel>_ALL_distances.csv
    <SeriesLabel>_summary_by_file.csv
    <SeriesLabel>_node_occupancy.csv
    ALL_pooled_by_series_summary.csv
    ALL_pooled_by_mRNA_summary.csv
    ALL_pooled_by_series_and_mRNA_summary.csv
```

### All-series outputs

```text
node_mrna_dist_BY_SERIES/ALL_SERIES_mRNA_to_nearest_node_distances.csv
node_mrna_dist_BY_SERIES/ALL_SERIES_summary_by_file.csv
node_mrna_dist_BY_SERIES/ALL_SERIES_node_occupancy.csv
node_mrna_dist_BY_SERIES/ALL_pooled_by_series_summary.csv
node_mrna_dist_BY_SERIES/ALL_pooled_by_mRNA_summary.csv
node_mrna_dist_BY_SERIES/ALL_pooled_by_series_and_mRNA_summary.csv
node_mrna_dist_BY_SERIES/ALL_rejected_coordinate_tables.csv
node_mrna_dist_BY_SERIES/ALL_unmatched_or_failed_files.csv
```

### Output interpretation

- `*_distances.csv`: one row per mRNA point, including mRNA coordinates, nearest-node coordinates, nearest node index, and distance in microns.
- `*_node_occupancy.csv`: one row per node, including how many mRNAs had that node as their nearest node.
- `*_summary.csv`: per-file summary of mean/median/min/max/std nearest-node distances plus counts within threshold distances.
- pooled summaries: aggregate across series and/or mRNA labels.

---

## 22. `mrna_to_nearest_node_distance_analysis.py` and older `mrna_to_nearest_node_distance_BY_SERIES.py`

These are older node-distance scripts. Keep them for provenance, but use the fixed BY_SERIES script above for new analyses.

### Older output folder names

```text
node_mrna_distance_output/
node_mrna_distance_output_BY_SERIES/
```

### Older output files

```text
ALL_mRNA_to_nearest_node_distances.csv
ALL_per_file_nearest_node_summary.csv
ALL_node_occupancy.csv
ALL_mRNA_nearest_node_distance_histogram.png
ALL_mRNAs_within_node_distance_thresholds.png
rejected_coordinate_tables.csv
unmatched_or_failed_files.csv
```

The fixed script renamed the main output folder to:

```text
node_mrna_dist_BY_SERIES
```

to reduce Windows path-length failures.

---

# PART F — Output map: what to look for in a completed folder

## 23. At the Series level

```text
Series XX/
├── MaxProjs.tif
├── RoiSet.zip
├── cells/
│   ├── single-cell TIFFs
│   ├── MitoGraph .txt skeleton files
│   ├── MitoGraph *_nodes.vtk files
│   ├── Hyperstacks_Grayscale/
│   └── mito_visualization_interactive/*.html
├── MS2/.../results/spots_extraction/
│   ├── converted_coordinates/cell_###_xyz_um.csv
│   ├── MS2_NN_distance_um.npy
│   ├── MS2_spot_counts.csv
│   ├── MS2_NN_distribution.png
│   ├── random_MS2_output/random_MS2_nn.npy
│   └── plots/MS2_real_vs_random.png
└── ATP2 or ATP3 or TIM50/.../results/spots_extraction/
    ├── converted_coordinates/cell_###_xyz_um.csv
    ├── <mRNA>_NN_distance_um.npy
    ├── <mRNA>_spot_counts.csv
    ├── <mRNA>_NN_distribution.png
    ├── random_<mRNA>_output/random_<mRNA>_nn.npy
    └── plots/<mRNA>_real_vs_random.png
```

## 24. At the condition parent level

```text
Condition parent folder/
├── NN_distance_recursive_run_summary.csv
├── random_mrna_recursive_run_summary_SHARED_SKELETON.csv
├── real_vs_random_recursive_comparison_summary.csv
├── pooled_real_vs_random_comparisons/
│   ├── pooled_input_file_pairing_summary.csv
│   ├── pooled_real_vs_random_summary.csv
│   ├── <mRNA>_POOLED_real_vs_random_REPLICATE_STATS.png
│   ├── <mRNA>_POOLED_real_NN_distance_um_RAW.npy
│   ├── <mRNA>_POOLED_random_NN_distance_um_RAW.npy
│   ├── <mRNA>_POOLED_real_NN_distance_um_DIFFRACTION_FLOORED.npy
│   ├── <mRNA>_POOLED_random_NN_distance_um_DIFFRACTION_FLOORED.npy
│   ├── <mRNA>_POOLED_real_vs_random_values.csv
│   └── <mRNA>_REPLICATE_LEVEL_real_vs_random_summary.csv
├── RNA_intensity_histograms/       optional, only if RNA intensity histogram script was run
└── node_mrna_dist_BY_SERIES/       optional, only if node-distance script was run
```


Optional RNA-intensity QC outputs, if `rna_intensity_histograms_recursive.py` was run:

```text
Condition or experiment parent folder/RNA_intensity_histograms/
├── RNA_intensity_histogram_analysis_workbook.xlsx
├── ALL_RNA_spot_intensities.xlsx
├── RNA_intensity_summary_by_series.xlsx
├── RNA_intensity_summary_by_series_and_channel.xlsx
├── RNA_intensity_summary_by_condition.xlsx
├── RNA_intensity_summary_by_condition_and_channel.xlsx
├── RNA_intensity_file_summary.xlsx
├── RNA_intensity_plot_index.xlsx
├── RNA_intensity_rejected_files.xlsx
├── RNA_intensity_RUN_SETTINGS.xlsx
└── plots/
    ├── parent/
    ├── by_series/
    ├── by_series_and_channel/
    ├── by_condition/
    └── by_condition_and_channel/
```

Optional mRNA-mRNA colocalization outputs, if `mrna_colocalization_by_cell.py` was run:

```text
Condition parent folder/mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/
├── ALL_cell_colocalization_summary.xlsx
├── ALL_<MRNA_A>_nearest_<MRNA_B>.xlsx
├── ALL_<MRNA_B>_nearest_<MRNA_A>.xlsx
├── ALL_<MRNA_A>_<MRNA_B>_pairs_within_<distance>.xlsx
├── ALL_<MRNA_A>_vs_<MRNA_B>_colocalization_workbook.xlsx
├── ALL_percent_colocalized_by_cell.png
├── ALL_<MRNA_A>_nearest_<MRNA_B>_distance_histogram.png
├── ALL_<MRNA_B>_nearest_<MRNA_A>_distance_histogram.png
└── per_series/<SeriesLabel>/cell_###/
```

Optional mRNA-to-mitochondrial-surface outputs, if `mito_rna_surface_GLOBAL_CALIBRATION.py` was run:

```text
Condition parent folder/mito_rna_surface_GLOBAL_CALIBRATION/
├── PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv
├── PARENT_spot_level_surface_distances_GLOBAL_CALIBRATION.csv
└── PARENT_channel_summary_GLOBAL_CALIBRATION.csv

Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/
├── GLOBAL_CALIBRATION_ACCEPTED_raw_spot_files.csv
├── GLOBAL_CALIBRATION_REJECTED_nonraw_files.csv
├── cell_###_<Channel>_GLOBAL_CALIBRATION_surface_distances.csv
└── cell_###_GLOBAL_CALIBRATION_surface_overlay.html
```

## 25. At the experiment parent level across strains/probes

```text
Experiment parent folder/
└── strain_probe_comparison/
    ├── strain_probe_folder_mapping.csv
    ├── strain_probe_series_level_metrics.csv
    ├── probe_set_mRNA_assignment_counts.csv
    ├── strain_probe_summary.csv
    ├── statistics_by_probe_across_strains_delta_median.csv
    ├── statistics_by_strain_across_probes_delta_median.csv
    ├── wt_vs_ko_statistics_delta_median.csv
    ├── strain_probe_delta_median_grouped_boxplot.png
    └── individual probe plots, if enabled
```

---

# PART G — How the scripts connect

## 26. Dependency chain

```text
smallFISH spots_extractions_*.csv
        ↓
rna_intensity_histograms_recursive.py  →  RNA_intensity_histograms/
        ↓
6.mrna_nn_distance_recursive.py
        ↓
converted_coordinates/cell_###_xyz_um.csv
<mRNA>_spot_counts.csv
<mRNA>_NN_distance_um.npy / .xlsx
        ↓                         ↓                         ↓
mrna_colocalization_by_cell.py     optional node-distance     7.random_mrna_null_recursive.py
        ↓                          scripts                   ↓
mrna_coloc_BY_SERIES/                                      random_<mRNA>_output/random_<mRNA>_nn.npy
                                                                ↓
                                                     8.real_vs_random_recursive_comparison.py
                                                                ↓
                                                     per-series plots and summary CSVs
                                                                ↓
                                                     9 or 10 pooled real-vs-random scripts
                                                                ↓
                                                     pooled_real_vs_random_comparisons/
                                                                ↓
                                                     11.compare_strains_probe_sets_boxplot_v6.py

MitoGraph *_mitosurface.vtk + raw smallFISH spots_extractions_*.csv
        ↓
mito_rna_surface_GLOBAL_CALIBRATION.py
        ↓
mito_rna_surface_GLOBAL_CALIBRATION/
```

## 27. Minimum files required for a completed real-vs-random condition

A condition is ready for pooled real-vs-random analysis when each mRNA has matched real and random files:

```text
<mRNA>_NN_distance_um.npy
random_<mRNA>_output/random_<mRNA>_nn.npy
```

and the parent folder contains the relevant `Series XX` folders.

A condition is ready for cross-strain/probe comparison when it contains:

```text
pooled_real_vs_random_comparisons/pooled_input_file_pairing_summary.csv
```

and the paired `.npy` files referenced by that CSV still exist or can be recovered by filename under the selected experiment parent folder.

---

# PART H — QC checklist


## 27A. After background subtraction

Check several representative TIFFs before continuing. Because `0.Recursive_Subtract_Background_OVERWRITE_ORIGINALS.ijm` overwrites the original images, compare against a backup raw image if possible.

Confirm that:

- the rolling-ball radius was `30` pixels unless there was a specific reason to change it,
- background haze is reduced,
- RNA puncta are still visible,
- mitochondrial signal is still continuous enough for MitoGraph,
- no channel appears over-subtracted or mostly black,
- the same TIFFs were not accidentally processed more than once.

---

## 28. Before MitoGraph

- Confirm cropped cells are centered and not clipped.
- Confirm mitochondrial channel has good signal-to-background.
- Confirm hyperstack dimensions are correct.
- Confirm channel extraction saved the correct channel names.

## 29. After MitoGraph

- Confirm `.txt` skeleton files exist in each `Series XX/cells` folder.
- Run `5.mito_vis_per_cells_folder_recursive.py`.
- Open a few `mito_visualization_interactive/*.html` files per series.
- Reject or rerun series with poor skeletonization.

## 30. After smallFISH

- Confirm every mRNA channel has a `results/spots_extraction` folder.
- Confirm cell-level coordinate tables exist.
- Confirm the detected spots visually match expected mRNA signal.



## 30A. After RNA intensity histogram script

Check for:

```text
RNA_intensity_histograms/RNA_intensity_histogram_analysis_workbook.xlsx
RNA_intensity_histograms/ALL_RNA_spot_intensities.xlsx
RNA_intensity_histograms/RNA_intensity_summary_by_series.xlsx
RNA_intensity_histograms/RNA_intensity_summary_by_condition.xlsx
RNA_intensity_histograms/RNA_intensity_plot_index.xlsx
RNA_intensity_histograms/plots/by_series/
RNA_intensity_histograms/plots/by_condition/
```

Open the histogram PNGs and check:

- intensity distributions are not unexpectedly shifted between replicate Series,
- one Series does not have an obvious low-intensity or saturated-intensity artifact,
- each RNA channel has a reasonable number of detected spots,
- `RNA_intensity_rejected_files.xlsx` does not show unexpected parsing failures,
- the run settings match the intended histogram bins and intensity-axis limits.

This step is especially useful before comparing spatial outputs across experimental conditions.

## 31. After real NN script

Check for:

```text
converted_coordinates/cell_###_xyz_um.csv
<mRNA>_NN_distance_um.npy
<mRNA>_spot_counts.csv
<mRNA>_NN_distribution.png
NN_distance_recursive_run_summary.csv
```

Open the run summary and check:

- `processed_files` is reasonable,
- `skipped_files` is not unexpectedly high,
- `pooled_nn_distances` is not zero unless expected.


## 32. After mRNA-mRNA colocalization script

Check for:

```text
mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/ALL_cell_colocalization_summary.xlsx
mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/ALL_<MRNA_A>_vs_<MRNA_B>_colocalization_workbook.xlsx
mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/ALL_unmatched_or_failed_cells.xlsx
```

Open the summary workbook and check:

- matched cell count is reasonable,
- `ALL_unmatched_or_failed_cells.xlsx` does not contain unexpected systematic missing cells,
- percent colocalized values are plausible,
- per-cell XY plots match visual expectations for a few representative cells.

## 33. After mitochondrial surface calibration script

For corrected datasets, prefer `mito_rna_surface_CORRECTED_XY_0p065_CROP250.py`. The older `GLOBAL_CALIBRATION` version is retained mainly for comparison with previous analyses.

Check for:

```text
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_spot_level_surface_distances_GLOBAL_CALIBRATION.csv
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_channel_summary_GLOBAL_CALIBRATION.csv
Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/cell_###_GLOBAL_CALIBRATION_surface_overlay.html
Series XX/cells/mito_rna_surface_CORRECTED_XY_0p065_CROP250/cell_###_CORRECTED_XY_0p065_CROP250_surface_overlay.html
```

Open several HTML overlays and check that RNA spots are properly aligned with the mitochondrial surface and skeleton. Also check the accepted/rejected raw spot file tables to confirm the script used raw smallFISH files, not generated downstream outputs.

## 34. After random null script

Check for:

```text
random_<mRNA>_output/random_<mRNA>_nn.npy
random_<mRNA>_output/random_<mRNA>_per_cell_summary.csv
random_mrna_recursive_run_summary_SHARED_SKELETON.csv
```

Open the per-cell summary and check:

- most cells have `Status = processed`,
- skeleton files were matched correctly,
- generated point counts match mRNA counts,
- missing skeletons are not systematic.

## 35. After pooled analysis

Check for:

```text
pooled_real_vs_random_comparisons/pooled_input_file_pairing_summary.csv
pooled_real_vs_random_comparisons/pooled_real_vs_random_summary.csv
```

Use `pooled_input_file_pairing_summary.csv` as the audit trail. It records which real and random files were paired, their series folder, status, and counts.

Use replicate-level metrics for inference. Do not rely only on pooled point-level p-values.

## 36. After cross-strain comparison

Check for:

```text
strain_probe_comparison/strain_probe_series_level_metrics.csv
strain_probe_comparison/strain_probe_summary.csv
strain_probe_comparison/wt_vs_ko_statistics_delta_median.csv
strain_probe_comparison/strain_probe_delta_median_grouped_boxplot.png
mrna_coloc_BY_SERIES/<comparison>/ALL_percent_colocalized_by_cell.png
Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/cell_###_GLOBAL_CALIBRATION_surface_overlay.html
Series XX/cells/mito_rna_surface_CORRECTED_XY_0p065_CROP250/cell_###_CORRECTED_XY_0p065_CROP250_surface_overlay.html
```

Each plotted dot should represent one series-level real-vs-random metric, not one individual mRNA spot.

---

# PART I — Common failure modes and fixes

## 37. No `spots_extraction` folders found

You selected the wrong parent folder or smallFISH outputs were not saved under the selected condition folder. Select the condition folder that contains all `Series XX` folders.



## 37A. RNA intensity histogram script finds no intensity values

Likely causes:

- selected parent folder does not contain raw `spots_extractions` files,
- smallFISH outputs were moved or renamed,
- files do not contain an `intensity` column,
- the CSV delimiter or output format changed,
- the files are inside a skipped generated-output folder,
- the selected folder is too narrow and does not include the Series or condition folders.

Check:

```text
RNA_intensity_histograms/RNA_intensity_file_summary.xlsx
RNA_intensity_histograms/RNA_intensity_rejected_files.xlsx
RNA_intensity_histograms/RNA_intensity_RUN_SETTINGS.xlsx
```

If the intensity column has a different name, set `INTENSITY_COLUMN` near the top of `rna_intensity_histograms_recursive.py` to the exact column name.


## 37B. Background subtraction is too weak, too strong, or accidentally repeated

Likely causes:

- rolling-ball radius was changed from the standard `30` pixels,
- the same TIFFs were processed more than once,
- raw images were not backed up before the overwrite macro was run,
- the selected parent folder contained already processed TIFFs,
- signal puncta are unusually broad and the rolling-ball radius is too small for that dataset.

Fix:

1. Go back to the backed-up raw TIFFs.
2. Rerun `0.Recursive_Subtract_Background_OVERWRITE_ORIGINALS.ijm` once.
3. Use `30` pixels as the standard rolling-ball radius.
4. Recheck a few RNA and mitochondrial channels before continuing.

Do not attempt to repair over-subtracted images by running later scripts. MitoGraph and smallFISH should be run from correctly background-subtracted TIFFs.

---

## 38. Real NN script makes no NN distances

Likely causes:

- all detected cells have fewer than two spots,
- coordinate columns were not detected,
- smallFISH output format changed,
- the wrong parent folder was selected.

Check:

```text
<mRNA>_NO_NN_DISTANCES.txt
NN_distance_recursive_run_summary.csv
```

## 39. Random null script reports missing shared skeletons

Likely causes:

- MitoGraph was not run for that series,
- skeleton `.txt` files are not in `Series XX/cells`,
- cell numbering between spot counts and skeleton files does not match,
- MitoGraph output files were moved outside the series folder.

Check:

```text
random_<mRNA>_per_cell_summary.csv
random_mrna_recursive_run_summary_SHARED_SKELETON.csv
```

## 40. Real-vs-random script cannot find random files

The comparison scripts expect:

```text
random_<mRNA>_output/random_<mRNA>_nn.npy
```

If you used the distance-range random script, the random file may be named:

```text
random_<mRNA>_nn_dist_<MIN>_to_<MAX>um.npy
```

Either use the standard random script for the primary workflow or update/rename the chosen random output with clear documentation.


## 41. mRNA colocalization script finds no matched cells

Likely causes:

- `6.mrna_nn_distance_recursive.py` was not run first,
- `converted_coordinates/cell_###_xyz_um.csv` files are missing,
- the two selected mRNAs do not both exist in the same Series folders,
- cell indices do not match between channels,
- `MRNA_A` or `MRNA_B` is misspelled.

Check:

```text
mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/ALL_input_coordinate_file_index.xlsx
mrna_coloc_BY_SERIES/<MRNA_A>_vs_<MRNA_B>_within_<distance>/ALL_unmatched_or_failed_cells.xlsx
```

## 42. Surface calibration script reports missing RNA files or missing mito surfaces

Likely causes:

- raw smallFISH files do not contain `spots_extractions` in the filename,
- channel folders are not named consistently with `RNA_CHANNELS`,
- MitoGraph did not create `*_mitosurface.vtk` files,
- `cells` folders were moved outside the Series folder,
- generated outputs were selected instead of the original condition parent folder.

Check:

```text
Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/GLOBAL_CALIBRATION_ACCEPTED_raw_spot_files.csv
Series XX/cells/mito_rna_surface_GLOBAL_CALIBRATION/GLOBAL_CALIBRATION_REJECTED_nonraw_files.csv
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv
```


## 42A. Corrected surface overlay is mirrored or shifted

If `mito_rna_surface_CORRECTED_XY_0p065_CROP250.py` still produces a visibly wrong overlay, troubleshoot in this order:

1. Confirm MitoGraph was run with the corrected microscope calibration, usually `-xy 0.065 -z 0.2`.
2. Confirm the cropped cell image size. The corrected script assumes `CROP_WIDTH_PIXELS = 250` and `CROP_HEIGHT_PIXELS = 250`.
3. If RNA points appear far too low or high, check `CROP_HEIGHT_PIXELS` first.
4. If RNA points are mirrored, test the transform setting:

```python
FIXED_TRANSFORM = "none"
FIXED_TRANSFORM = "flip_y"
FIXED_TRANSFORM = "flip_x"
FIXED_TRANSFORM = "flip_xy"
```

5. If the overlay is close but uniformly shifted, adjust only:

```python
GLOBAL_DX_UM
GLOBAL_DY_UM
GLOBAL_DZ_UM
```

Only change one category at a time: calibration, crop size, transform, or offset. After each change, rerun one representative cell and inspect `cell_###_CORRECTED_XY_0p065_CROP250_surface_overlay.html` before processing the full dataset.

---

## 43. Windows path-length errors

Use the short-path versions of scripts where available. The updated scripts save short cell labels such as:

```text
cell_000
cell_001
```

instead of repeating entire ND2 filenames in output paths.

The fixed node-distance script also changed the output folder from:

```text
node_mrna_distance_output_BY_SERIES
```

to:

```text
node_mrna_dist_BY_SERIES
```

## 44. Node-distance script matches only `cell_002` or mismatches `_000_nodes.vtk`

Use:

```text
mrna_to_nearest_node_distance_BY_SERIES_FIXED.py
```

The fixed script matches `cell_###_xyz_um.csv` to `*_###_nodes.vtk` using the final node-file suffix rather than unrelated numbers in the ND2 filename.

## 45. Derived surface-distance CSVs are treated as input mRNA files

Use the fixed node-distance script. It restricts input discovery to:

```text
converted_coordinates/cell_*_xyz_um.csv
```

and excludes:

```text
mito_rna_surface_GLOBAL_CALIBRATION/*_surface_distances.csv
```
---

## 46. Using different probes other than MS2, ATP2, ATP3, or TIM50

Future datasets may use probes that are not part of the current default set, such as a new nuclear-encoded mitochondrial mRNA, a different endogenous transcript, or a different reporter/probe combination. The workflow can still be used, but the new probe name must be added anywhere the Python scripts use a fixed list of recognized mRNA/RNA channel names.

### Why this matters

Several scripts infer the RNA identity from folder names and filenames. For example, the script may decide that a file belongs to `MS2`, `ATP2`, `ATP3`, or `TIM50` by searching the path for one of those strings. If a new probe is named `COX1`, `COX2`, `VAR1`, or another transcript name, the script may label it as `unknown`, skip it, or fail to match it to the correct cell unless that new name is added to the recognized-probe list.

### Use consistent folder names first

Before editing scripts, make sure the folder naming is consistent across all Series folders. Use the exact same spelling and capitalization wherever possible.

Good example:

```text
Series 4/COX1/.../results/spots_extraction/spots_extractions_..._000.csv
Series 5/COX1/.../results/spots_extraction/spots_extractions_..._000.csv
```

Avoid mixing names such as:

```text
Cox1
COX-1
COXI
cox1_probe
```

unless the scripts are explicitly updated to recognize each version.

### Where to add a new probe name

Open each Python script in a text editor and search for lists like these:

```python
KNOWN_MRNA_NAMES = ["MS2", "ATP2", "ATP3", "TIM50"]
KNOWN_RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "TIM50", "Tim50"]
RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "Tim50"]
```

Add the new probe name to the list. For example, if the new probe is `COX1`, change the list to:

```python
KNOWN_MRNA_NAMES = ["MS2", "ATP2", "ATP3", "TIM50", "COX1"]
```

or:

```python
RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "Tim50", "COX1"]
```

Do this in every script that needs to recognize the new probe.

### Scripts that commonly need this update

Check these scripts first:

- `6.mrna_nn_distance_recursive.py`: search for `KNOWN_MRNA_NAMES` or the supported mRNA-name list. This lets the script label the new probe correctly and write `<probe>_NN_distance_um.npy`, `<probe>_NN_distance_um.xlsx`, and `<probe>_spot_counts.csv/.xlsx`.

- `7.random_mrna_null_recursive.py`: search for `KNOWN_MRNA_NAMES` or the mRNA-name list. This lets the random null model find the new probe's spot-count file and write `random_<probe>_output/`.

- `7.1 random_mrna_distance_range.py`: search for `KNOWN_MRNA_NAMES` or the mRNA-name list. This allows distance-range randomization for the new probe.

- `8.real_vs_random_recursive_comparison.py`: search for the mRNA-name list or filename-matching section. This lets real and random `.npy` files for the new probe pair correctly.

- `9.pooled_real_vs_random_across_series_REPLICATE_LEVEL_STATS.py`: search for the mRNA-name list or parsing function. This allows pooled real-vs-random analysis to include the new probe.

- `10. pooled_real_vs_random_ Thesis Graphs.py`: search for the mRNA-name list or parsing function. This allows the new probe to appear in thesis-style pooled figures.

- `11.compare_strains_probe_sets_boxplot_v6.py`: search for the probe-set assignment section. This is needed if the new probe should appear as its own category in the strain/probe comparison plots.

- `mrna_colocalization_by_cell.py`: search for `KNOWN_MRNA_NAMES`, `MRNA_A`, and `MRNA_B`. This is needed to compare the new probe against MS2, ATP2, ATP3, TIM50, or another probe.

- `mito_rna_surface_GLOBAL_CALIBRATION.py`: search for `RNA_CHANNELS`. This is needed so the new probe is included in mRNA-to-mitochondrial-surface proximity analysis.

- `rna_intensity_histograms_recursive.py`: search for `KNOWN_RNA_CHANNELS`. This is needed so RNA intensity histograms are grouped under the correct probe name instead of `unknown_RNA`.

- `mrna_to_nearest_node_distance_BY_SERIES_FIXED.py`: search for `KNOWN_MRNA_NAMES`. This is needed so mRNA-to-node distance outputs are labeled by the correct probe.

### Example: adding COX1 to the workflow

If the new probe is `COX1`, update the relevant script settings like this:

```python
KNOWN_MRNA_NAMES = ["MS2", "ATP2", "ATP3", "TIM50", "COX1"]
```

For the mitochondrial surface script:

```python
RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "Tim50", "COX1"]
```

For the RNA intensity histogram script:

```python
KNOWN_RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "TIM50", "Tim50", "COX1"]
```

For the mRNA-mRNA colocalization script, also change the pair being compared:

```python
MRNA_A = "MS2"
MRNA_B = "COX1"
COLOCALIZATION_DISTANCE_UM = 0.25
```

### Check output names after adding a probe

After updating the scripts, the new probe should appear in output filenames and summary tables. For example, `COX1` should generate files such as:

```text
COX1_NN_distance_um.npy
COX1_NN_distance_um.xlsx
COX1_spot_counts.csv
COX1_spot_counts.xlsx
random_COX1_output/random_COX1_nn.npy
RNA_intensity_summary_by_condition_and_channel.xlsx  # with COX1 rows
```

If the new probe appears as `unknown`, `unknown_RNA`, or is missing from summaries, the script did not recognize the probe name. Re-open the script and search again for the recognized-name list or path-parsing function.

### Update cross-strain/probe comparisons carefully

The cross-strain comparison script may have custom logic that groups probes into biological probe sets. For example, ATP6 and ATP8 may be collapsed into `ATP6/8`, and MS2 rows may be interpreted as the ATP6/8 reporter. If a new probe should be plotted as its own category, add it to the probe-assignment logic and confirm it appears correctly in:

```text
strain_probe_comparison/probe_set_mRNA_assignment_counts.csv
strain_probe_comparison/strain_probe_series_level_metrics.csv
strain_probe_comparison/strain_probe_summary.csv
```

Do not assume the new probe will automatically appear in the final strain/probe figure just because earlier scripts processed it. The final comparison script may still need an explicit probe-set mapping.

### Quick troubleshooting checklist

If a new probe does not work, check the following in order:

1. The raw smallFISH folder is named consistently with the probe name.
2. The `spots_extractions_*.csv` or `.xlsx` files exist inside the expected `results/spots_extraction/` folder.
3. The new probe name was added to every relevant `KNOWN_MRNA_NAMES`, `KNOWN_RNA_CHANNELS`, or `RNA_CHANNELS` list.
4. The script was saved after editing.
5. The selected parent folder contains the Series folders, not just one output subfolder.
6. Audit files such as `*_input_coordinate_file_index.xlsx`, `RNA_intensity_file_summary.xlsx`, or accepted/rejected raw spot file tables show the new probe under the correct name.

When in doubt, search the script for `MS2`, `ATP2`, `ATP3`, `TIM50`, `Tim50`, `KNOWN_MRNA_NAMES`, `KNOWN_RNA_CHANNELS`, and `RNA_CHANNELS`. Any section that explicitly lists old probe names may need the new probe added.

---

# PART J — Final handoff notes

For long-term storage, each completed analysis folder should contain:

1. the raw or preprocessed image folders,
2. `Series XX` folders,
3. MitoGraph outputs,
4. smallFISH outputs,
5. script-generated summaries,
6. pooled comparison outputs,
7. cross-condition comparison outputs, if applicable,
8. this procedure document.

The most important audit-trail CSV files are:

```text
NN_distance_recursive_run_summary.csv
random_mrna_recursive_run_summary_SHARED_SKELETON.csv
real_vs_random_recursive_comparison_summary.csv
pooled_real_vs_random_comparisons/pooled_input_file_pairing_summary.csv
pooled_real_vs_random_comparisons/pooled_real_vs_random_summary.csv
strain_probe_comparison/strain_probe_series_level_metrics.csv
strain_probe_comparison/strain_probe_summary.csv
RNA_intensity_histograms/RNA_intensity_histogram_analysis_workbook.xlsx
RNA_intensity_histograms/RNA_intensity_summary_by_series.xlsx
RNA_intensity_histograms/RNA_intensity_summary_by_condition.xlsx
mrna_coloc_BY_SERIES/<comparison>/ALL_cell_colocalization_summary.xlsx
mrna_coloc_BY_SERIES/<comparison>/ALL_<MRNA_A>_vs_<MRNA_B>_colocalization_workbook.xlsx
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv
mito_rna_surface_GLOBAL_CALIBRATION/PARENT_channel_summary_GLOBAL_CALIBRATION.csv
```

The most important final figure files are:

```text
pooled_real_vs_random_comparisons/<mRNA>_POOLED_real_vs_random_REPLICATE_STATS*.png
strain_probe_comparison/strain_probe_delta_median_grouped_boxplot.png
RNA_intensity_histograms/plots/by_series/*.png
RNA_intensity_histograms/plots/by_condition/*.png
```

The most important reusable numeric arrays are:

```text
<mRNA>_NN_distance_um.npy
random_<mRNA>_output/random_<mRNA>_nn.npy
<mRNA>_POOLED_real_NN_distance_um_RAW.npy
<mRNA>_POOLED_random_NN_distance_um_RAW.npy
<mRNA>_POOLED_real_NN_distance_um_DIFFRACTION_FLOORED.npy
<mRNA>_POOLED_random_NN_distance_um_DIFFRACTION_FLOORED.npy
```

If a person opens this folder later, these files should let them reconstruct the entire path from raw images to final figures and strain/probe statistical summaries.
