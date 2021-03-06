%% Purpose
% This is a wrapper script for creating the voxel-wise abnormality maps of patients, in
% the units of standard deviation
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs 
% CAMinputs.patientCode=Unique ID of the patients
% CAMinputs.pathOfMetricImages= path of the normative database
% CAMinputs.pathOfPatientParamImg=path of the CMRGlc images
% CAMinputs.pathOfT1MRImg=Path of the t1-mprage mr image - original resolution.
% CAMinputs.pathOfDartelTemplate=path of the dartel template, which was generated by rundartel sub-routine.
% CAMinputs.pathOfAbnormalityMaps=path of the abnormality maps - where it should be stored?
% CAMinputs.pathOfHS=path of hammersmith atlas

%%Inputs from the user: 

% !! Please fill in the required fields accurately !!

CAMinputs.patientCode='';
CAMinputs.pathOfMetricImages= '';
CAMinputs.pathOfPatientParamImg='';
CAMinputs.pathOfT1MRImg='';
CAMinputs.pathOfDartelTemplate='';
CAMinputs.pathOfAbnormalityMaps='';
CAMinputs.pathOfHS='';

%% Run the program
createAbnormalityMap(CAMinputs);