%% Purpose
% This is a wrapper script for running the 'performMRbasedMOCO.m' program, 
% which is responsible for alignment of MR masks with the reconstructed PET
% frames
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs to be entered by the user                  
%mrMoCoInputs.pathOfDicomPET = Path of reconstructed dicom pet
%mrMoCoInputs.pathOfMRnavigators = Path of MR navigators
%mrMoCoInputs.pathOfDICOMMRmask=Path of segmented MR masks
%mrMoCoInputs.pathToStoreMoCoAnalysis=path to store the MR moco analysis
%mrMoCoInputs.pathOfT1mr=path of T1 MR
%mrMoCoInputs.pathOfDixonInPhase=path of Dixon in-phase
%mrMoCoInputs.pathOfTOFMRA=path of TOF-MRA
%mrMoCoInputs.pathOfNiftiPET = Path of SPM converted nifty pet.

%% Limitations
% Only "DICOM" images are allowed - NIFTI support not included at this
% point in time, but will be added in due course.

%% Program start
% Copy your physical path of the 3D time-of-flight MR angiography DICOM
% images.
mrMoCoInputs.pathOfDicomPET = [pathPrefix,volunteerID,'/Processed data/Reconstructed PET/MoCo_Recons_PET_CT'];
mrMoCoInputs.pathOfMRnavigators =[pathPrefix,volunteerID,'/Processed data/MR navigators'];
mrMoCoInputs.pathOfDICOMMRmask=[pathPrefix,volunteerID,'/Processed data/MR scaled petrous'];
mrMoCoInputs.pathToStoreMoCoAnalysis=[pathPrefix,volunteerID,'/Processed data'];
mrMoCoInputs.pathOfT1mr=[pathPrefix,volunteerID,'/Processed data/T1 MR'];
mrMoCoInputs.pathOfDixonInPhase=[pathPrefix,volunteerID,'/Processed data/Attenuation correction maps/Static AC/DIXON/HEAD_FDG_PET_MRAC_LM_60MIN_IN_0009'];
mrMoCoInputs.pathOfTOFMRA=[pathPrefix,volunteerID,'/Processed data/MR TOF'];
mrMoCoInputs.pathOfNiftiPET =[pathPrefix,volunteerID,'/Processed data/Reconstructed PET/Nifti_MoCo_Recons_PET_CT'];
performMRbasedMOCO(mrMoCoInputs);

 