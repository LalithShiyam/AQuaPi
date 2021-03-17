%% Purpose
% This is a wrapper script for running the 'createDartelNormDB.m' program, 
% which is responsible for creating normative database.
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs 
% * DartelDBinputs.pathOfParamImg -> Path of the folder where the parametric images are stored.
% * DartelDBinputs.pathOfCoregT1  -> Path of the folder where the coregistered T1 MR images are stored.

%% Limitations
% All the images should be in nifti or analyse format.

%% Inputs  to be entered by the user
 DartelDBinputs.pathOfParamImg='';  % please enter the respective physical paths
 DartelDBinputs.pathOfCoregT1='';
 
 %% Program start
 createDartelNormDB(DartelDBinputs);
 