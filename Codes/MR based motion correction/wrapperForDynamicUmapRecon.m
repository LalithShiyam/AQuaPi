%% Purpose
% This is a wrapper script for running the 'coreDynamicUmapRecon.m' program, 
% which is responsible reconstructing the PET list-mode after the alignment
% of u-maps with the PET emission data.
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs to be entered by the user                  
%  Inputs: 
%       DURinputs.subjectID= ID of the subject.
%       DURinputs.pathOfStaticACmaps= path of the static u-maps (can be
%       more than 1)
%       DURinputs.pathOfReconParamFile= path of the reconstruction
%       parameter file
%       DURinputs.pathOfDicomMRnavigators= path where the DICOM MR
%       navigators are stored (EPI sequences need to be in 'split' mode)
%       DURinputs.pathToListModePET= path pointing towards the PET
%       list-mode data.
%       DURinputs.pathToCopyReconPETdata = path to copy the reconstructed
%       pet data.
%       
%  Outputs: 
%       Dynamic reconstructions of the PET data with motion induced AC
%       maps.

%% Limitations
% Only "SPLIT" EPI images are allowed - Mosaic images support not included at this
% point in time, but will be added in due course.

%% Program start

% Fill in the required details

DURinputs.subjectID= '';
DURinputs.pathOfStaticACmaps= '';
DURinputs.pathOfReconParamFile= '';
DURinputs.pathOfDicomMRnavigators= '';
DURinputs.pathToListModePET= '';
DURinputs.pathToCopyReconPETdata = '';
disp('Performing dynamic u-map PET LM reconstruction using JS-recon...')
coreDynUmapReconstruction(DURinputs)

%%


