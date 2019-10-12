%% Purpose
% This is a wrapper script for running the 'convertCTtoMRumap.m' program,
% which converts a low-dose CT image to Siemens compatible CT-AC map.
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs to be entered by the user                  
%         umapGen.pathOfLowDoseCT - path to the low dose CT dicom series
%         umapGen.pathOfDixonSeries - path to the Dixon series (F, IN, OPP, W, U-map)
%         umapGen.pathToStoreUmap - where to store the u-maps.

%% Limitations
% Only "DICOM" images are allowed - NIFTI support not included at this
% point in time, but will be added in due course.

%% Program start
% Copy your physical path of the 3D time-of-flight MR angiography DICOM
% images.

umapGen.pathOfLowDoseCT=''; % To be filled in by the user!
umapGen.pathOfDixonSeries = ''; % To be filled in by the user!
umapGen.pathToStoreUmap = ''; % To be filled in by the user!

convertCTtoMRumap(umapGen); % Running the program.