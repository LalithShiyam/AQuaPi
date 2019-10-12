%% Purpose
% This is a wrapper script for running the 'coreSegmentation.m' program, 
% which is responsible for automatic petrous and cervical
% segmentation.
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs to be entered by the user                  
% coreSegInputs.pathOfAngio - Physical path of the 3D time-of-flight MR
% angiography DICOM images.

%% Limitations
% Only "DICOM" images are allowed - NIFTI support not included at this
% point in time, but will be added in due course.

%% Program start
% Copy your physical path of the 3D time-of-flight MR angiography DICOM
% images.

coreSegInputs.pathOfAngio = ''; % To be filled in by the user!
coreSegInputs.patientCode = ''; % To be filled in by the user!
%% Hard-coded variables.

% Error messages

errorMsg{1}='This program needs a path to the "DICOM" series!';
errorMsg{2}='Raise an issue in github...';

% path to store the data.

cd(coreSegInputs.pathOfAngio); 
cd .. 
coreSegInputs.path2StoreSeg=pwd; % path where the output will be stored.


%% Preliminary checks are being done here, to see if the program can be run.
fileFormat=checkFileFormat(coreSegInputs.pathOfAngio); % check if the folder has dicom images.
switch fileFormat 
    case 'DICOM'
        disp(['DICOM images found in ',coreSegInputs.pathOfAngio,'...']);
        disp('Applying segmentation algorithm on the dataset...');
        coreSegmentation(coreSegInputs); % Running the coreSegmentation algorithm
    otherwise
        error(errorMsg{1});
        error(errorMsg{2});
end

%%


