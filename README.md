# Sample code for: Pallidal alpha activity is an electrophysiological biomarker of symptom severity in obsessive-compulsive disorder

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit)

## Introduction
This repository contains the code for the preprint [Vissani et al 2026](https://doi.org/10.21203/rs.3.rs-9316385/v1) "Pallidal alpha activity is an electrophysiological biomarker of symptom severity in obsessive-compulsive disorder". The full raw dataset is is available upon request.

Abstract of the paper:
>Clinical efficacy of deep brain stimulation (DBS) for obsessive–compulsive disorder (OCD) remains limited by the absence of objective neural biomarkers to guide targeting, programming, and longitudinal optimization of therapy. Using sensing-enabled DBS devices, we analyzed intracranial local field potentials recorded longitudinally from thirteen patients with treatment-refractory OCD implanted in the ventral capsule/ventral striatum (VC/VS). We identify the periodic component of alpha-band activity (7–14 Hz) as a robust electrophysiological correlate of OCD symptom severity. At the initial clinical evaluation, periodic alpha power explained inter-individual variance in symptom severity, and within individuals it tracked longitudinal symptom fluctuations during chronic treatment. Spatial mapping localized the source of this periodic alpha activity to the anterior external globus pallidus (GPe), aligning with stimulation contacts selected as therapeutically optimal during clinical monopolar review. Together, these findings identify anterior GPe alpha oscillations as a state-dependent and spatially specific biomarker of OCD severity, with direct implications for electrophysiology-informed DBS programming and the development of adaptive neuromodulation strategies.

This repository can be downloaded by entering the following commands:
```
git clone https://github.com/Brain-Modulation-Lab/Vissani2026_alphaOCDBiomarker.git
cd Vissani2026_alphaOCDBiomarker
```
## Minimal Dataset 

The minimum dataset `Vissani2026_alphaOCDBiomarker/demos` to run `ShowSpectralParameterization.m`, `ShowMonopolarEstimation.m` and `ShowSpatialAnalysis.m` contains:
* `sub-OP1002_ses-postop01_lead-dbs-mni-coords_hemi-left.csv` contains the Lead=DBS MNI coordinates of the DBS lead contacts
* `sub-OP1002_ses-postop01_percept-bss_hemi-left.mat` contains the electrophysiological neural signals recorded during the first recording session

## MATLAB Analysis

* SPECTRAL PARAMETERIZATION: the script `ShowSpectralParameterization.m` is designed to illustrate the computation of the absolute, relative and periodic power spectra
* PSEUDO-MONOPOLAR ESTIMATION: the script `ShowMonopolarEstimation.m` is designed to illustrate the computation of the estimation of the pseudo-monopolar contribution of power at each contact from the observation of bipolar signals.
* SPATIAL ANALYSIS: the script `ShowSpatialAnalysis.m` runs a demo that illustrates the computation of spatial focality in simulated data. The weighted spatial inertia is a concept from physics that quantifies the spread of a mass distribution. In this context, the power value served as the weight, and the inertia was calculated as:


```math
I = sum(w_i * d_i^2) / sum(w_i)
```

where:
- $w_i$ = power at contact i
- $d_i$ = distance from weighted centroid

The denominator ensures that has units of $mm^2$ and enables scale-invariant comparison across power definitions. Moreover, the square root of can be used to define the radius of gyration, representing the effective radius of an equivalent sphere that has the same moment of inertia as the observed spatial distribution.

## R Analysis
* The script `LongitudinalModellingBiomarker.R` runs the longitudinal modeling of OCD symptom severity based on electrophysiological features recorded over multiple recording sessions.

If you encounter any problems, please report them as issues in the repository or send an [email](mailto:mvissani@mgh.harvard.edu).
This repository has been tested successfully in MATLAB versions 2022a and 2023a on MacOS and Windows.


## External dependencies

The code depends on these repositories:

* [fieldtrip](https://www.fieldtriptoolbox.org/): toolbox to analyze electrophysiological data
* [bml](https://github.com/Brain-Modulation-Lab/bml): fieldtrip wrapper developed by the [BrainModulation Lab](https://www.brainmodulationlab.org/).

You need to manually download and include (only the main folder!) them in your MATLAB dependencies.
After that just run these commands in MATLAB to manage dependencies:
```
bml_defaults
ft_defaults
```

The external folder in the repo contains other libraries:
* [RainCloud](https://github.com/RainCloudPlots/RainCloudPlots): toolbox to visualize data distributions.
* [PERMUTOOLS](https://github.com/mickcrosse/PERMUTOOLS/tree/master): toolbox to implement permutation-based statistics.


## Contributors
* [Matteo Vissani](mailto:mvissani@mgh.harvard.edu)

>Citation: [INSERT DOI HERE]

## Full Dataset request
The full raw dataset is available upon request.

## Funding
This work was funded by a Michael Jenike Young Investigator Award of the International Obsessive Compulsive Disorder Foundation (IOCDF)

# License
This project is covered under the **MIT License**.
