%% Purpose
% This is a wrapper script for running the 'coreSegmentation.m' program, 
% which is responsible for automatic petrous and cervical
% segmentation.
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs                   
% coreSegInputs.pathOfAngio - Physical path of the 3D time-of-flight MR
% angiography DICOM images.

%% Limitations
% Only "DICOM" images are allowed - NIFTI support not included at this
% point in time, but will be added in due course.

%% Program start

% Copy your physical path of the 3D time-of-flight MR angiography DICOM
% images.

coreSegInputs.pathOfAngio = ''; % To be filled in by the user!

% Running the coreSegmentation algorithm

coreSegmentation(coreSegInputs);
