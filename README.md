# Cardiac-T2-Bull-s-Eye-Segmentation-and-AHA-Regional-Quantification

MATLAB workflow for interactive MRI analysis using ROI-based segmentation, AHA-style bull’s-eye mapping, QA overlays, and segment-wise export. Showcasing image processing, manual annotation, geometric region mapping, visualization, debugging, and reproducible quantitative analysis pipeline design.

## Core technical components

This project demonstrates practical implementation of:

- **MATLAB-based image analysis**
- **Interactive ROI drawing and manual annotation tools**
- **Cardiac slice and phase selection workflows**
- **Geometric region mapping and sector-based segmentation**
- **AHA-style myocardial segment labeling**
- **Quantitative extraction of regional imaging values**
- **Quality-control overlay generation**
- **Figure export and structured result saving**
- **Iterative debugging and version-controlled workflow refinement**
- **Reproducible analysis pipeline design**

## Project overview

This project was developed to generate a structured **bull’s-eye representation of myocardial T2\*** values from cardiac MRI data stored in `.mat` format. The workflow begins with time-resolved cardiac imaging data and allows the user to interactively select representative **basal, mid, and apical** short-axis slices, along with a chosen **cardiac phase** from the acquired sequence.

Using magnitude images as the primary anatomical reference, the user manually defines the **endocardial** and **epicardial** borders of the myocardium and marks a reference insertion point used to orient the myocardial segmentation. From these inputs, the script constructs a myocardial mask, assigns regional segment labels, and summarizes the T2\* values within each anatomical region.

The final output is a **bull’s-eye diagram**, a standard cardiac visualization format used to display regional myocardial measurements in a circular arrangement that reflects the left ventricle from base to apex. In addition to the bull’s-eye, the workflow also produces **QA overlay images**, segment-wise result tables, and MATLAB output files containing intermediate and final analysis products.

## Cardiac background

In cardiac imaging, the **left ventricle** is commonly divided into standardized anatomical regions so that regional tissue measurements can be compared consistently across slices, phases, and subjects. A bull’s-eye plot is a compact way to visualize these regional values by arranging them in concentric rings:

- the **outer ring** represents the **basal** myocardium
- the **middle ring** represents the **mid-ventricular** myocardium
- the **inner ring** represents the **apical** myocardium

In this project, the myocardium is organized according to an **AHA-style regional segmentation framework**, where the basal and mid levels are divided into six regions each and the apical level is divided into four. This allows segment-wise T2\* values to be extracted and visualized in a format that is easier to interpret than slice-wise raw maps alone.

T2\* mapping is useful because it provides a quantitative measure derived from MRI signal behavior, and regional analysis can help identify spatial variation across the myocardium rather than relying on only a global average.

## Workflow summary

The typical processing flow is:

1. Load a `.mat` file containing:
   - `T2_starmaps_all`
   - `iField`

2. Display magnitude images for anatomical guidance.

3. Select representative:
   - **basal**
   - **mid**
   - **apical**
   slices.

4. Select a cardiac phase from the available time frames.

5. Draw:
   - endocardial contour
   - epicardial contour

6. Mark the reference insertion point used for segment orientation.

7. Build the myocardial mask and assign regional segment labels.

8. Compute segment-wise T2\* statistics.

9. Export:
   - QA overlay image
   - bull’s-eye figure
   - segment-wise CSV table
   - result `.mat` files

## Input data

The workflow is designed for `.mat` files containing cardiac MRI-derived data, specifically variables such as:

- `T2_starmaps_all`
- `iField`

These inputs are used to generate:
- T2\* maps for quantitative analysis
- magnitude images for anatomical ROI definition and QA

## Output files

Depending on the version used, the pipeline can generate:

- **QA overlay images** showing ROI boundaries, segment labels, and reference points
- **bull’s-eye diagrams** summarizing regional myocardial values
- **CSV / Excel-compatible tables** containing segment-wise values
- **MATLAB `.mat` result files** storing masks, segment maps, and quantitative outputs
- compact analysis files containing:
  - `metricMaps`
  - `segLabelMaps`

## Version development

This project was developed iteratively to improve both the analysis logic and the interpretability of the outputs.

### Initial working version
The first working version established an end-to-end pipeline from raw `.mat` input to a basic bull’s-eye result using one representative basal, mid, and apical slice and one selected cardiac phase.

### Workflow expansion
Later versions extended the design toward a more flexible workflow by supporting broader slice handling, cleaner processing structure, and improved export logic.

### AHA-oriented refinement
Subsequent versions focused on aligning the segmentation more closely with standardized myocardial regional structure, improving the interpretation of the overlay image, refining the segment orientation logic, and stabilizing the full export pipeline.

### QA improvements
Later refinements added:
- clearer QA overlays
- better segment labeling
- visible segment boundaries
- more readable figure layouts
- cleaner visual presentation for supervisor review and validation

## Why this project matters

This project is not only a visualization tool but also a reproducible image-analysis workflow that combines:

- manual anatomical input
- quantitative extraction
- structured regional mapping
- interpretable visual outputs

It is useful as a research workflow for testing how regional myocardial values can be extracted from raw mapping data and translated into a standardized summary view that is easier to review, compare, and communicate.

## Skills reflected in this project

Without relying on black-box automation, this work highlights hands-on implementation in:

- computational image analysis
- research-focused MATLAB development
- user-guided annotation workflows
- segmentation logic design
- quantitative regional measurement
- scientific visualization
- exportable reporting
- iterative debugging and refinement of analysis software

## Repository contents

A typical repository structure may include:

```text
scripts/
docs/
results_example/
README.md
