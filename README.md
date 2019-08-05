![AQuaPi-Logo](AQuaPi-Logo.png)

# AQuaPi Toolbox

AQuaPi Suite: Absolute Quantification Suite, is a MATLAB based computational toolbox, for non-invasively measuring cerebral metabolic rates of glucose using the synergistic data from a fully-integrated PET/MR.The framework comprises of an image-derived input function component and a quantification component, that work together to produce cerebral metabolic rates of glucose maps.

# Abilities

- Petrous segmentation from a 3D time-of-flight MR angiography
- MR driven motion correction (motion navigators needed): alignment of attenuation maps possible.
- Iterative partial volume correction aware of spatiotemporal variations of the target and background activities
- Voxelwise Patlak analysis
- Spatial normalisation (SPM12 DARTEL)
- Z-score calculation

# Requirements

- MATLAB version > or = R2016a 
- SPM 12 
- Siemens Biograph mMR and e7 reconstruction tools
  - 3D time-of-flight MR angiography sequence (FOV dependent)
  - 60-min PET list-mode acquisition
  - 3D MPRAGE T1 MR sequence
  - 3D EPI MR navigators  (need to be in 'split' mode using Siemens SyngoVia)
  - Attenuation map: preferrably CT or pseudo-CT

# Literature

- Shiyam Sundar LK, Muzik O, Beyer T. (2018). Towards quantitative [18F]FDG-PET/MRI of the brain: Automated MR-driven calculation of an image-derived input function for the non-invasive determination of cerebral glucose metabolic rates. Journal of Cerebral Blood Flow & Metabolism. https://doi.org/10.1177/0271678X18776820
- Shiyam Sundar LK, Muzik O, Beyer T. (2019). The promise of fully-integrated PET/MR imaging: Non-invasive clinical quantification of cerebral glucose metabolism. The Journal of Nuclear Medicine. https://doi.org/10.2967/jnumed.119.229567
