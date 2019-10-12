%% Purpose
% This is a wrapper script for running the 'performPatlak.m' program, 
% which is responsible for voxel-wise patlak analylsis
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs to be entered by the user                  

%plasmaGlucose= plasma glucose levels in units mmol/L
%patlakInputs.cutOff=  time in seconds.
% patlakInputs.voxelWise= True or False (has regional patlak capabilities
% as well)
%patlakInputs.file2Delete='';
%patlakInputs.subjectID='';
%patlakInputs.pathOfNiftiPET=path of reconstructed nifty pet
%patlakInputs.pathOfTransformedMRnavigators=path of the transformed MR
%navigators (pet-space navigators)
%patlakInputs.pathOfMRcorrespondenceFile=path of mr correspondence file
%patlakInputs.plasmaGlucose=path of plasma glucose levels;
%patlakInputs.lumpedConstant=0.65;
%patlakInputs.pathToStore=where to store the images.
%patlakInputs.plasmaInputFunction=path where the plasma input function is
%present.
%patlakInputs.BrainMask=path of the brain masks.

%% Program start
% Copy your physical path of the 3D time-of-flight MR angiography DICOM
% images.
plasmaGlucose='' % units mmol/L
patlakInputs.cutOff= '' % time in seconds.
patlakInputs.voxelWise=1;
patlakInputs.file2Delete='';
patlakInputs.subjectID='';
patlakInputs.pathOfNiftiPET='';
patlakInputs.pathOfTransformedMRnavigators='';
patlakInputs.pathOfMRcorrespondenceFile='';
patlakInputs.plasmaGlucose='';
patlakInputs.lumpedConstant=0.65;
patlakInputs.pathToStore='';
patlakInputs.plasmaInputFunction='';
patlakInputs.BrainMask='';
performPatlak(patlakInputs);

