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
%mrMoCoInputs.pathOfNiftiPET = Path of reconstructed SPM converted nifty pet using
%"constructDynPETseries.m"

%% Limitations
% Only "DICOM" images are allowed - NIFTI support not included at this
% point in time, but will be added in due course.

%% Program start
% Copy your physical path of the images
pathPrefix='';
volunteerID='';
mrMoCoInputs.pathOfDicomPET = '';
mrMoCoInputs.pathOfMRnavigators ='';
mrMoCoInputs.pathOfDICOMMRmask='';
mrMoCoInputs.pathToStoreMoCoAnalysis='';
mrMoCoInputs.pathOfT1mr='';
mrMoCoInputs.pathOfDixonInPhase='';
mrMoCoInputs.pathOfTOFMRA='';
mrMoCoInputs.pathOfNiftiPET ='';
coreMRbasedMOCO(mrMoCoInputs); % Running the core MR based moco program.

 