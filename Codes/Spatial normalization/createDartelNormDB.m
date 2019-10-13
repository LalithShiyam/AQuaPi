%% Normative database creation

%% Author information
% 
%  Lalith Kumar Shiyam Sundar, M.Sc.
%  Quantitative Imaging and Medical Physics, Medical University of Vienna
%  Date: 27.11.2017
% 
% 
% * Wrapper program or main program: createDartelNormDB.m
% * Sub-routines: segmentForDartel_job.m | runDARTEL_job.m |
% normalizeToMNI_job.m | NormDBmetric.m

%% Purpose of the program
% This is a function which creates a DARTEL (spm 12) based normative
% database for the CMRGlc images (parametric images). The program performs
% the following steps in a streamlined manner:

%%
% * Unified segmentation: splitting the coregistered T1-MR images to DARTEL
% compatible gray matter and white matter. Also the inverse (MNI to patient
% space) and forward (patient to MNI space) transformations are generated 
% in this step (takes a lot of time).
% * Run DARTEL: in this step, the individual gray/white matter are warped
% with each other iteratively and the deformation fields are generated
% here along with the final template (template 6.nii), which is our
% reference template.
% * Normalize to MNI: here the deformation fields generated from above for 
% individual subjects will be used to warp the patient space parametric
% images, coregistered T1 images to the template space (which is in MNI)

%% Inputs for the program 
% * DartelDBinputs.pathOfParamImg -> Path of the folder where the parametric images are stored.
% * DartelDBinputs.pathOfCoregT1  -> Path of the folder where the coregistered T1 MR images are stored.

%% Outputs for the program 
% * Mean Image (MeanImg.nii) -> The average image based on the input
% parametric images
% * Standard deviation image (StdDevImg.nii) -> The standard deviation
% image bsaed on the input parametric images
% 

%% Unified segmentation
% Here, the coregistered T1 MR images of the subjects will be segmented
% into gray matter (has a prefix 'c1' for native space and 'rc1' for 
% DARTEL), white matter (has a prefix 'c2' for native space and 'rc2' for
% DARTEL). Others have a default setting, do not worry. Also the forward
% (patient space to MNI) and the inverse (MNI to patient space) are stored
% here in '.nii' format (NIFTI). 
%% 
% * Sub-routine script: segmentForDartel_job.m
% * Inputs: Coregistered T1 MR images filenames
% * Outputs: Deformation fields (forward: prefix 'y', inverse: prefix 'iy')
%            and DARTEL files (prefix 'r*');

function []=createDartelNormDB(DartelDBinputs)

cd(DartelDBinputs.pathOfCoregT1) % go to the folder containing the coregistered T1 MR images.
CoregMRfiles=dir('*Coreg*.nii'); % read the file names in the directory
for lp=1:length(CoregMRfiles) % prepare the input for the sub-routine
    CoregMRstrings{lp,:}=[pwd,filesep,CoregMRfiles(lp).name];
end
SFDinputs.CoregMRfilepath=CoregMRstrings; % input strings pushed to the sub-routine's input structure.

% Load tissue probability maps by searching the system
tpmFolderPath=what('tpm');
for lp=1:6 % the value is always 6 in tpm
    SFDinputs.tpm{lp,:}=[tpmFolderPath.path,filesep,'TPM.nii',',',num2str(lp)];
end
segmentForDartel_job(SFDinputs) % sub-routine function which runs the segmentation
disp('Gray and white matter segments for DARTEL created...')

%% Run DARTEL
% Here the gray and white matter of the control subjects are registered
% between themselves in an iterative manner, the final output is a template
% which is in the MNI space. This template can be considered as an
% average brain generated from your personal control group.
%%
% 
% * Sub-routine script: runDARTEL_job.m
% * Inputs: DARTEL gray (prefix 'rc1') and white matter (prefix 'rc2')
% * Outputs: template_6 -> which is the average brain generated based on
%            your control group brain and deformation fields for individual
%            subjects (prefix 'u*');
%
% 
cd(DartelDBinputs.pathOfCoregT1)% go to the folder containing the coregistered T1 MR images.
DARTELgm=dir('rc1*.nii'); % read in the DARTEL gray matter inputs
DARTELwm=dir('rc2*.nii'); % read in the DARTEL white matter inputs
for lp=1:length(DARTELgm)
    tempDARTELgmfiles{lp,:}=[pwd,filesep,DARTELgm(lp).name]; % load the path for DARTEL compatible gray matter
end
for lp=1:length(DARTELwm)
    tempDARTELwmfiles{lp,:}=[pwd,filesep,DARTELwm(lp).name]; % load the path for DARTEL compatible white matter
end
runDARTELinputs.GrayWhiteCell={tempDARTELgmfiles,tempDARTELwmfiles}; % load the paths in a single cell and flush it to the input for the sub-routine
%runDARTEL_job(runDARTELinputs); % create the template
disp('DARTEL template and deformation fields created...')

%% Normalize to MNI space
% After the creation of the control group based template (template_6.nii),
% all the control subject's images (parametric images + coregistered T1)
% need to be warped to the MNI space, this can be done by using the
% deformation fields (prefix 'u*') and the final template.
%%
% 
% * Sub-routine script: normalizeToMNI_job.m
% * Inputs: DARTEL generated flow fields (prefix 'u*'), average brain 
%           template generated from DARTEL (template_6.nii), coregistered
%           T1 MR images and parametric images of the control subjects in
%           their native space.
% * Outputs: Normalized parametric images and coregistered T1 MR images in
%           MNI space
% 
cd(DartelDBinputs.pathOfCoregT1)% go to the folder containing the coregistered T1 MR images.
TemplateOfInterest=dir('Template_6.nii'); % always load the last template, this is average sharp template
normToMNIinputs.pathToTemplate{1}=[pwd,filesep,TemplateOfInterest.name]; % save the path of the template to the input structure 
flowFieldsFiles=dir('u*.nii'); %find the file names which correspond to the flow fields
for lp=1:length(flowFieldsFiles)
    normToMNIinputs.pathToFlowFields{lp,:}=[pwd,filesep,flowFieldsFiles(lp).name]; % flush the flow field paths to the input structure
end
CoregMRfiles=dir('Coreg*.img'); %load the coregistered T1 MR images
for lp=1:length(CoregMRfiles)
    tempMRfiles{lp,:}=[pwd,filesep,CoregMRfiles(lp).name];
end
cd(DartelDBinputs.pathOfParamImg) % go to the folder where the parametric images in native space are present.
ParamImgFiles=dir('*.nii'); % load all the images with the extension '*.img')
for lp=1:length(ParamImgFiles)
    tempParamFiles{lp,:}=[pwd,filesep,ParamImgFiles(lp).name];
end
normToMNIinputs.pathToImgToNormalize={tempMRfiles,tempParamFiles}; % load the paths in a single cell and flush it to the input for the sub-routine
normalizeToMNI_job(normToMNIinputs); % feed the input to the sub-routine and watch it run.
disp('Normative database created...')

%% Metric images creation
% Here the mean and the standard deviation images are created based on the
% normalized control group parametric images. Note that the parametric
% images which have been normalized are in MNI space.
%%
% * Sub-routine script: NormDBmetric.m
% * Inputs: Parametric images normalized to the MNI space
% * Outputs: Mean image, standard deviation image based on the control
%            group of normalized parametric images (MNI)
cd(DartelDBinputs.pathOfParamImg) % go to the folder where the parametric images in native space are present.
cd ..
mkdir('Metric Images')
DartelDBinputs.pathOfMetricImg=[pwd,filesep,'Metric Images']; % create the path for the metric images and flush it to the input parameter of the global function.
cd(DartelDBinputs.pathOfParamImg) % go to the folder where the parametric images in native space are present.
NormFileNames=dir('sw*.nii'); % load the MNI normalized parametric images (prefix 'sw')
for lp=1:length(NormFileNames)
    fNames{lp,:}=NormFileNames(lp).name; % input for the file names
end
MetricImg=NormDBmetric(fNames,DartelDBinputs.pathOfMetricImg); % Mean and standard deviation images are created and written to the directory.
disp('Mean and standard deviation images created!')
end
