%% Creation of abnormality maps and inverse mapping (DARTEL -> Individual space)
% This function basically creates the patient abnormality maps in units of
% standard deviation (Z-maps) by comparing the patient specific parametric
% image with the control database. Note that all the images are normalized
% to DARTEL space. To make it clinically usable, we inverse map the
% abnormality map from the DARTEl space to its corresponding individual
% space.
% 
%  Author
%  Lalith Kumar Shiyam Sundar, M.Sc.
%  Quantitative Imaging and Medical Physics, Medical University of Vienna
%  30.11.2017
%  Routine name: createAbnormalityMap.m
%  Sub-routines:

%% Inputs
% The inputs for this particular program is going to be the path where the
% metric images are stored, the path of patient parametric images, the path
% of the T1 MR images and the path where the template_6.nii is
% stored.
%patientList={'UT'};
%patientList={'BM','CS','DV','HA','KN','KS','PD','RS','SE','UT' };
function []=createAbnormalityMap(CAMinputs)
leftPutamenIndex=38;rightPutamenIndex=39;
    %% Coregister the patient parametric image with its T1 MR image.
    % This is mandatory, make sure you do this. DO NOT SKIP IT!

    cd(CAMinputs.pathOfT1MRImg); % go to the path containing the T1 MR image and make sure that the T1 Image is the only image in the folder.
    convertDicomtoNii(cd,cd);
    fName=dir('T1*.img') % Look for analyze images.
    if isempty(fName)
        fName=dir('T1*.nii'); % if empty look for nifty images.
    end
    CoregInputs.SourceImgPath=[pwd,filesep,fName.name]; % T1 MR image passed to the coregistration job module (moving image)
    CoregInputs.Interp=1; 
    CoregInputs.Prefix='Coreg_';
    cd(CAMinputs.pathOfPatientParamImg) % go to the folder containing the parametric image, make sure this particular folder has just one parametric image.
    fName=[]; % recycling
    fName=dir('*.img') % Look for analyze images.
    if isempty(fName)
        fName=dir('*.nii'); % if empty look for nifty images.
    end
    CoregInputs.RefImgPath=[pwd,filesep,fName.name]; % pass the parametric image to the coregistration job module (stationary image)
    cd(CAMinputs.pathOfT1MRImg); % go to the path containing the T1 MR image, just to put the MR results in a single place.
    CoregInputs.MaskImgPath={''}; % nothing to coregister in tandem.
   Coregistration_job(CoregInputs) % run the coregistratoin module

    %% Run segment (split the gray and white matter of the patient T1 MR)
    % search for SPM
    WhereIsSPM12=what('spm12');
    for lp=1:6
    TPMstringsForSegment{lp,:}=[WhereIsSPM12.path,filesep,'tpm',filesep,'TPM.nii',',',num2str(lp)];
    end
    SFDinputs.tpm=TPMstringsForSegment;
    cd(CAMinputs.pathOfT1MRImg); % making sure we are still in the T1 MR folder
    fName=[];
    fName=dir('Coreg*.nii');
    SFDinputs.CoregMRfilepath={[pwd,filesep,fName.name]}; % flush the Coregistered MR file path to the segment module.
   segmentForDartel_job(SFDinputs); % run the segmentation

    %% Run DARTEL (existing template)
    % This module is responsible for bringing the segmented gray and white
    % matter to the DARTEL template space.
    cd(CAMinputs.pathOfT1MRImg)% go to the folder containing the coregistered T1 MR images.
    DARTELgm=dir('rc1*.nii'); % read in the DARTEL gray matter inputs
    DARTELwm=dir('rc2*.nii'); % read in the DARTEL white matter inputs
    for lp=1:length(DARTELgm)
        tempDARTELgmfiles{lp,:}=[pwd,filesep,DARTELgm(lp).name]; % load the path for DARTEL compatible gray matter
    end
    for lp=1:length(DARTELwm)
        tempDARTELwmfiles{lp,:}=[pwd,filesep,DARTELwm(lp).name]; % load the path for DARTEL compatible white matter
    end
    REDinputs.GrayWhiteCell={tempDARTELgmfiles,tempDARTELwmfiles}; % load the paths in a single cell and flush it to the input for the sub-routine
    cd(CAMinputs.pathOfDartelTemplate) % go to the path containing the templates
    tmplFnames=dir('Template*.nii') % load the template file names;
    tmplFnames=tmplFnames(2:end); % remove template_0.nii
    for lp=1:length(tmplFnames)
        REDinputs.pathOfTemplateSeries{lp}=[pwd,filesep,tmplFnames(lp).name]
    end
   runExistingDartel_job(REDinputs) % run the DARTEL tool for existing modules.
    disp('DARTEL template and deformation fields created for patients...')

    %% Normalize the patient parametric image to MNI space with all the deformations
    cd(CAMinputs.pathOfT1MRImg); % go to the folder containing the coregistered  T1 MR images.
    flowFieldsFiles=dir('u*.nii'); %find the file names which correspond to the flow fields
    cd(CAMinputs.pathOfDartelTemplate);% go to the folder path containing DARTEL template
    TemplateOfInterest=dir('Template_6.nii'); % always load the last template, this is average sharp template
    normToMNIinputs.pathToTemplate{1}=[pwd,filesep,TemplateOfInterest.name]; % save the path of the template to the input structure 
    cd(CAMinputs.pathOfT1MRImg); % go to the folder containing the coregistered  T1 MR images.
    for lp=1:length(flowFieldsFiles)
        normToMNIinputs.pathToFlowFields{lp,:}=[pwd,filesep,flowFieldsFiles(lp).name]; % flush the flow field paths to the input structure
    end
    CoregMRfiles=dir('Coreg*.img'); %load the coregistered T1 MR images
    if isempty(CoregMRfiles)
        CoregMRfiles=dir('Coreg*.nii'); % look for nifty if analyze is not available in the folder.
    end
    for lp=1:length(CoregMRfiles)
        tempMRfiles{lp,:}=[pwd,filesep,CoregMRfiles(lp).name];
    end
    grayWhiteDartelFiles=dir('rc*.nii');
    if isempty(grayWhiteDartelFiles)
        grayWhiteDartelFiles=dir('rc*.img');
    end
    for lp=1:length(grayWhiteDartelFiles)
        tempGrayWhite{lp,:}=[pwd,filesep,grayWhiteDartelFiles(lp).name];
    end
    cd(CAMinputs.pathOfPatientParamImg) % go to the folder where the parametric images in native space are present.
    ParamImgFiles=dir('*.img'); % load all the images with the extension '*.img')
    if isempty(ParamImgFiles)
        ParamImgFiles=dir('*.nii'); % look for nifty, if analyze not available in the folder.
    end
    for lp=1:length(ParamImgFiles)
        tempParamFiles{lp,:}=[pwd,filesep,ParamImgFiles(lp).name];
    end
    normToMNIinputs.pathToImgToNormalize={tempMRfiles,tempParamFiles,{tempGrayWhite{1}},{tempGrayWhite{2}}}; % load the paths in a single cell and flush it to the input for the sub-routine
    normalizeToMNI_job(normToMNIinputs); % feed the input to the sub-routine and watch it run.
    disp('Patient parametric image normalized in DARTEL Space...') % output will be 'sw' prefix files.

    %% Load the metric images : mean and standard deviation

    cd(CAMinputs.pathOfMetricImages) % go to the folder containing the metric images
    stdImg=spm_read_vols(spm_vol('StdDevImg.nii'));
    meanImg=spm_read_vols(spm_vol('MeanImg.nii'));
    meanImgMask=meanImg>=0;
    meanImg=meanImg.*meanImgMask;
    
    
   
    
    

    %% Generation of the abnormality maps

    dummyHeader=spm_vol('MeanImg.nii');
    dummyHeader.fname=['Abnormality Map patient ',CAMinputs.patientCode,'.nii'];
    dummyHeader.pinfo=[1;0;0];
    cd(CAMinputs.pathOfT1MRImg);
    grayProbMask=dir('swrc1*.nii');% read the normalized gray matter
    whiteProbMask=dir('swrc2*.nii');% read the normalized white matter
    grayProbMask=spm_read_vols(spm_vol(grayProbMask.name));
    whiteProbMask=spm_read_vols(spm_vol(whiteProbMask.name));
    BrainMask=grayProbMask+whiteProbMask;
    BrainMask=BrainMask>0.10;
    BrainMask=imfill(BrainMask,'holes'); % create a brain mask from the patient itself
    cd(CAMinputs.pathOfPatientParamImg); % go to the path of the parametric image
    ImgOfInterestFile=dir('sw*.nii'); % get the resliced parametric image of the patient in DARTEL standard space.
    ImgOfInterest=spm_read_vols(spm_vol(ImgOfInterestFile.name));
    ImgOfInterestMask=ImgOfInterest>=0;
    ImgOfInterest=ImgOfInterestMask.*ImgOfInterest;
    AbnormalityMap=round((ImgOfInterest-meanImg)./stdImg).*BrainMask;
    HyperActiveAreas=((AbnormalityMap>0).*AbnormalityMap).*BrainMask;
    HypoActiveAreas=abs((AbnormalityMap<0).*AbnormalityMap).*BrainMask;
    AbnormalityMap(isnan(AbnormalityMap))=0;
    AbnormalityMap(isinf(AbnormalityMap))=0;
    HypoHdr=dummyHeader;
    cd(CAMinputs.pathOfAbnormalityMaps);
    HypoHdr.fname=['HypoActiveArea_',CAMinputs.patientCode,'.nii'];
    HyperHdr=dummyHeader;
    HyperHdr.fname=['HyperActiveArea_',CAMinputs.patientCode,'.nii'];
    spm_write_vol(dummyHeader,AbnormalityMap);
    spm_write_vol(HyperHdr,HyperActiveAreas);
    spm_write_vol(HypoHdr,HypoActiveAreas);
    disp('Abnormality maps for patients created...');
   
    %% Crude patch: Normalizing to putamen - correcting for physiological state:
   	
    cd(CAMinputs.pathOfPatientParamImg); % go to the path of the parametric image
    ImgOfInterestFile=dir('sw*.nii'); % get the resliced parametric image of the patient in DARTEL standard space.
    ImgOfInterest=spm_read_vols(spm_vol(ImgOfInterestFile.name));
    ImgOfInterestMask=ImgOfInterest>=0;
    cd(CAMinputs.pathOfHS);
    HSfile=dir('*.hdr');
    V=spm_vol(HSfile.name);
    voxsiz = [1.5 1.5 1.5]; % new voxel size {mm}
    %John gems to reslice objects to defined voxel size.
    for i=1:numel(V)
       bb        = spm_get_bbox(V(i));
       VV(1:2)   = V(i);
       VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
       VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
       VV(1).dim = VV(1).dim(1:3);
       spm_reslice(VV,struct('mean',false,'which',1,'interp',0)); % 1 for linear
    end
    reslicedHS=dir('r*.hdr')
    hammersmithAtlas=spm_read_vols(spm_vol(reslicedHS.name));
    newHdr=spm_vol(reslicedHS.name);
    leftPutamenAtlas=hammersmithAtlas==leftPutamenIndex;
    leftPutamen=ImgOfInterest(leftPutamenAtlas);
    rightPutamenAtlas=hammersmithAtlas==rightPutamenIndex;
    rightPutamen=ImgOfInterest(rightPutamenAtlas);
    patientPutamen=median([leftPutamen;rightPutamen]);
    normalPutamen=median([meanImg(rightPutamenAtlas);meanImg(leftPutamenAtlas)]);
    scalingFactor=normalPutamen./patientPutamen;
    ImgOfInterest=(ImgOfInterestMask.*ImgOfInterest).*scalingFactor;
    normAbnormalityMap=round((ImgOfInterest-meanImg)./stdImg).*BrainMask;
    HyperActiveAreas=((normAbnormalityMap>0).*normAbnormalityMap).*BrainMask;
    HypoActiveAreas=abs((normAbnormalityMap<0).*normAbnormalityMap).*BrainMask;
    normAbnormalityMap(isnan(normAbnormalityMap))=0;
    normAbnormalityMap(isinf(normAbnormalityMap))=0;
    HypoHdr=dummyHeader;
    cd(CAMinputs.pathOfAbnormalityMaps);
    HypoHdr.fname=['NormToPutamen_HypoActiveArea_',CAMinputs.patientCode,'.nii'];
    HyperHdr=dummyHeader;
    HyperHdr.fname=['NormToPutamen_HyperActiveArea_',CAMinputs.patientCode,'.nii'];
    normParamImgHdr=HyperHdr;
    normParamImgHdr.fname=['NormParamImg_',CAMinputs.patientCode,'.nii'];
    spm_write_vol(normParamImgHdr,ImgOfInterest);
    spm_write_vol(dummyHeader,normAbnormalityMap);
    spm_write_vol(HyperHdr,HyperActiveAreas);
    spm_write_vol(HypoHdr,HypoActiveAreas);
    
    disp('Putamen scaled Abnormality maps for patients created...');
    
     
    %% Invert the abnormality maps from the DARTEL/MNI space to individual space
    cd(CAMinputs.pathOfT1MRImg); % go to the MR path=
    deformFile=dir('u_r*.nii');
    InvDartelInputs.pathOfDartelFlowField=[pwd,filesep,deformFile.name];
    InvDartelInputs.pathOfTemplate=[CAMinputs.pathOfDartelTemplate,filesep,'Template_6.nii'];
    CoregT1File=dir('Coreg*.nii');
    InvDartelInputs.pathOfBaseImg=[pwd,filesep,CoregT1File.name];
    cd(CAMinputs.pathOfMetricImages) % go to the folder containing the metric images
    copyfile('*.nii',CAMinputs.pathOfAbnormalityMaps);
    cd(CAMinputs.pathOfAbnormalityMaps);
    fNames=dir('*.nii');
    for lp=1:length(fNames)
        InvDartelInputs.pathOfFilesToInvert{lp,:}=[pwd,filesep,fNames(lp).name];
    end
    cd(CAMinputs.pathOfAbnormalityMaps)
    mkdir('Inverse warped')
    cd([pwd,filesep,'Inverse Warped']);
    InvDartelInputs.pathToStoreInverseImages=cd;
    InverseDartel_job(InvDartelInputs)
    disp('Inverse warping from DARTEL to individual space done!');
    
end
