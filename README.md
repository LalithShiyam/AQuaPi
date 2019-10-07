![AQuaPi-Logo](AQuaPi-Logo.png)

# AQuaPi Toolbox

AQuaPi: Absolute Quantification Pipeline, is a MATLAB based computational toolbox, for non-invasively measuring cerebral metabolic rates of glucose using the synergistic data from a fully-integrated PET/MR.The framework comprises of an image-derived input function component and a quantification component, that work together to produce cerebral metabolic rates of glucose maps.

Table of Contents
=================

  * [Abilities](#Abilities)
  * [Requirements](#Requirements)
  * [Literature](#Literature)
  * [Installation](#Installation)
  * [Usage](#Usage)
    * [Automated segmentation](#Automated-segmentation)
    * [MR-driven motion compensation](#MR-driven-motion-compensation)
    * [Iterative regional partial volume correction](#Iterative-regional-partial-volume-correction)
    * [Voxel-wise Patlak analysis](#Voxel-wise-Patlak-analysis)
    * [Spatial normalisation](#Spatial-normalisation)
    * [Z-Score calculation](#Z-Score-calculation)
  * [Contributing](#Contributing)
  * [Credits](#Credits)
  * [License](#License)
  


# Abilities

- Petrous/cervical segmentation from a 3D time-of-flight MR angiography dataset
- MR driven motion correction (motion navigators needed): alignment of attenuation maps possible.
- Iterative regional partial volume correction aware of the activity variations of the target and background.
- Voxelwise Patlak analysis
- Spatial normalisation (SPM12 DARTEL)
- Z-score calculation

# Requirements

- Software requirements
  - MATLAB version > or = R2016a 
  - SPM 12 
  - Siemens e7 reconstruction tools + JS-recon
  - Powershell ISE > 4.0
- Input data: Refer to the individual pages for the suggested settings and corresponding limitations
  - 3D time-of-flight MR angiography
  - 60-min PET list-mode acquisition
  - 3D MPRAGE T1 MR sequence
  - 3D EPI MR navigators
  - Attenuation map: preferrably CT or pseudo-CT

# Literature

- Shiyam Sundar LK, Muzik O, Beyer T. (2018). Towards quantitative [18F]FDG-PET/MRI of the brain: Automated MR-driven calculation of an image-derived input function for the non-invasive determination of cerebral glucose metabolic rates. Journal of Cerebral Blood Flow & Metabolism. https://doi.org/10.1177/0271678X18776820
- Shiyam Sundar LK, Muzik O, Beyer T. (2019). The promise of fully-integrated PET/MR imaging: Non-invasive clinical quantification of cerebral glucose metabolism. The Journal of Nuclear Medicine. https://doi.org/10.2967/jnumed.119.229567

# Installation

# Usage

  # Automated segmentation

