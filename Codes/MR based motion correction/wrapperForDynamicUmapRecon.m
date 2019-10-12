
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
% Only CT/Dixon/UTE/PseudoCT/Resolute Rotation is possible at this point, support to other AC
% maps will be implemented in due course.



%% Program start

% Fill in the required details

DURinputs.subjectID= 'HC002';
DURinputs.pathOfStaticACmaps= '/Users/lalith/Documents/HC002/AC maps/CTAC'; 
DURinputs.pathOfReconParamFile= '/Users/lalith/Documents/HC002/AC maps/Reconstruction parameter file';
DURinputs.pathOfDicomMRnavigators= '/Users/lalith/Documents/HC002/AC maps/MR navigators';
DURinputs.pathToListModePET= '/Users/lalith/Documents/HC002/AC maps/LM PET';
DURinputs.pathToCopyReconPETdata = '/Users/lalith/Documents/HC002/AC maps/Recon PET';
disp('Performing dynamic u-map PET LM reconstruction using JS-recon...')
coreDynUmapReconstruction(DURinputs)

%%


