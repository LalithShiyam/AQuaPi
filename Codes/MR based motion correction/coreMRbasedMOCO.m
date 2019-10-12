%-------------------------------------------------------------------------%
% This is a function which alings the MR time of angiography mask to the
% PET frames, based on the motion vector fields obtained from the MR
% navigators (which are approximated to be in correspondence with their
% respective PET frames based on least temporal difference)
% 
% Author: Lalith Kumar Shiyam Sundar, M.Sc.
% Date: Feb 05, 2018
% Quantitative Imaging and Medical Physics
% 
% Inputs: 
%         [1] mrMoCoInputs.pathOfDicomPET = path of the dicom PET data (one
%           frame is enough)
%         [2] mrMoCoInputs.pathOfMRnavigators= path of the split dicom MR navigators
%         [3] mrMoCoInputs.pathOfDICOMMRmask = path of the DICOM MR mask, which is
%           in alignment with the first MR navigator.
%         [4] mrMoCoInputs.pathToStoreMoCoAnalysis = path to store the
%           entire analysis. 
%         [5] mrMoCoInputs.pathOfT1mr= path of the T1 MR.  
%         [6] mrMoCoInputs.pathOfDixonInPhase = path of dixon in-phase u-map
%         [7] mrMoCoInputs.pathOfTOFMRA = path of 3d-tof mr angiography.
%         [8] mrMoCoInputs.pathOfNiftiPET= path of the Nifti PET images (for performing PET based MoCo as well : incase mr based MoCo fails)  
% Output:
%         - "Dynamic MR masks" folder, which contains the same number of MR
%         rotated MR masks as the MR navigators 
% 
%
% Pseudoalgorithm: 
% -----------------
%        - Read PET dicom images.
%        - Read split MR navigators and resample them to PET dicom space
%        (the split  MR navigators can be fused with pet images due to voxel correspondence)
%        - 
%
% Usage: performMRbasedMoCo(mrMoCoInputs);
%
%
%
%
%
%                   
%-------------------------------------------------------------------------%
%                           Program Start 
%-------------------------------------------------------------------------%

function []=coreMRbasedMOCO(mrMoCoInputs)

%% Create a folder called MoCo
cd(mrMoCoInputs.pathToStoreMoCoAnalysis);
WhereAmI=cd; % finding out where i am.
mkdir('MoCo') % creating a folder called MoCo. 
pathOfWorkingFolder=[pwd,filesep,'MoCo'];
%% Go to the folder containing the DICOM PET data.

cd(mrMoCoInputs.pathOfDicomPET); % go to the path containing the dicom PET. (chances are that the folder is split)

% find if there are any subfolders
subFolders=dir;
subFolders=subFolders([subFolders(:).isdir]==1);
subFolders=subFolders(~ismember({subFolders(:).name},{'.','..'}));
if isempty(subFolders)
   PET=readImages(cd,'align');
else
   temp=[pwd,filesep,subFolders(1).name];
   cd(temp);
   PET=readImages(cd,'align');
end
disp('DICOM PET files read to workspace...')

%% Go to the folder containing the split MR navigators in DICOM format.

cd(mrMoCoInputs.pathOfMRnavigators)
navFolders=dir('*SPLIT*');
navFolders=navFolders(arrayfun(@(x) x.name(1),navFolders) ~= '.');

% sort the MR navigator folders using "natsort"

for lp=1:length(navFolders) 
    foldersToSort{lp,:}=navFolders(lp).name;
end

sortedFolders=natsort(foldersToSort);
wb=timebar('Transforming images from MR space to PET-space ...');

parfor lp=1:length(foldersToSort)
    cd([mrMoCoInputs.pathOfMRnavigators,filesep,sortedFolders{lp}])
    MRnavigator=readImages(cd,'align');
    outMRnav=transformVoxelSpace(MRnavigator,PET);
    cd(pathOfWorkingFolder); %% change it to the MoCo path
    writeDicomImage(outMRnav,['Transform_MRnav_',num2str(lp),'_To_PET',]);
    timebar(wb,lp/length(foldersToSort));
end
close(wb);
disp('Transforming MR navigators to PET space...');



%% Convert the transformed MR navigators and MR mask to NIFTI format and store them in a seperate subfolder inside MoCo folder

cd(pathOfWorkingFolder)
pathOfNifti=[pwd,filesep,'Nifti working folder'];
subFolders=dir;
subFolders=subFolders(arrayfun(@(x) x.name(1), subFolders) ~= '.'); % read volunteer folders
mkdir('Nifti working folder')
parfor lp=1:length(subFolders)
    temp=[pwd,filesep,subFolders(lp).name];
    cd(temp);
    convertDicomtoNii(cd,cd);
    oldNiftiFile=dir('*nii');
    movefile(oldNiftiFile.name,[pathOfNifti,filesep,subFolders(lp).name,'.nii'])
    cd ..
end
disp('Converted all the DICOM images to nifti format...')
%% Go the folder containing dixon in-phase series 

cd(mrMoCoInputs.pathOfDixonInPhase);
inPhaseDixon=readImages(cd);
transInPhaseDixon=transformVoxelSpace(inPhaseDixon,PET);
cd(pathOfWorkingFolder)
writeDicomImage(transInPhaseDixon,'Transform_Dixon_In_Phase_PET');
pathOfTransDixon=[pathOfWorkingFolder,filesep,'Transform_Dixon_In_Phase_PET'];
convertDicomtoNii(pathOfTransDixon,pathOfTransDixon);
oldFile=dir('*nii');
movefile(oldFile.name,[pathOfNifti,filesep,'Transform_Dixon_In_Phase_PET.nii']);
stationaryImg=[pathOfNifti,filesep,'Transform_Dixon_In_Phase_PET.nii'];

% load mr tof 

cd(mrMoCoInputs.pathOfTOFMRA)
convertDicomtoNii(cd,cd);oldFile=dir('*nii');
movefile(oldFile.name,[pathOfNifti,filesep,'TOF_MRA.nii']);
movingImg=[pathOfNifti,filesep,'TOF_MRA.nii'];

% load MR mask

cd(mrMoCoInputs.pathOfDICOMMRmask);
convertDicomtoNii(cd,cd);
oldFile=dir('*nii');
movefile(oldFile.name,[pathOfNifti,filesep,'MR_mask.nii']);
pathOfMRmask={[pathOfNifti,filesep,'MR_mask.nii']};


% coregistration between tof-mra (moving) and dixon in-phase (Stationary),
% petrous mask (mr mask - other image)

CoregInputs.Interp=0; % nearest neighbor interpolation, as its a binary mask
CoregInputs.Prefix='Coreg_';
CoregInputs.RefImgPath=stationaryImg; 
CoregInputs.SourceImgPath=movingImg;
CoregInputs.MaskImgPath=pathOfMRmask;
Coregistration_job(CoregInputs);

cd(pathOfNifti)
mrMaskFile=dir('Coreg*MR*mask*');
movefile(mrMaskFile.name,[pathOfNifti,filesep,'MR_mask_1.nii']);
%delete MR_mask.nii

%% Go to the folder containing T1 MR volume

cd(mrMoCoInputs.pathOfT1mr)
convertDicomtoNii(cd,cd);
oldNiftiFile=dir('*nii');
movefile(oldNiftiFile.name,[pathOfNifti,filesep,'T1_MR','.nii']);

%% segmentation of the t1-MR to obtain the brain mask. 

cd(pathOfNifti)
T1MR=dir('T1_MR*');
segmentBrainInputs.pathOfT1MR=[pathOfNifti,filesep,T1MR.name];
WhereIsSPM12=what('spm12');
for lp=1:6
TPMstringsForSegment{lp,:}=[WhereIsSPM12.path,filesep,'tpm',filesep,'TPM.nii',',',num2str(lp)];
end
segmentBrainInputs.TPMstringPath=TPMstringsForSegment;
segmentBrain(segmentBrainInputs)

%% Create a brain mask and feed it as another mask to the rotation program (below).

cd(pathOfNifti)
brainFiles=dir('c*T1*MR*')
for lp=1:length(brainFiles)
    sortingBrainFiles{lp,:}=brainFiles(lp).name;
end
sortedBrainFiles=natsort(sortingBrainFiles);
brainMask=zeros(size(spm_read_vols(spm_vol(sortedBrainFiles{1,:}))));
for lp=1:2 % including just grey and white matter. no csf.
    brainMask=brainMask+spm_read_vols(spm_vol(sortedBrainFiles{lp,:}));
end
brainMask=brainMask>0.45; % converting the brain mask into a logical.

dummyHeader=spm_vol(sortedBrainFiles{1,:});
dummyHeader.pinfo=[1;0;0];
dummyHeader.fname='Brain_Mask_1.nii';
spm_write_vol(dummyHeader,brainMask);
dummyHeader.fname='Brain_mask.nii';
spm_write_vol(dummyHeader,brainMask);

%% Perform coregistration between T1 MR (source - moving) and Coreg-TOF MR-angiography (reference - fixed) 

CoregInputs.Interp=0; % Nearest neigbour interpolation
cd(pathOfNifti);
movingImg=dir('T1_MR*nii');
CoregInputs.SourceImgPath=[pathOfNifti,filesep,movingImg.name];
cd(pathOfNifti);
CoregInputs.MaskImgPath={[pathOfNifti,filesep,'Brain_Mask_1.nii']}; % no mask at this point
fixedImg=dir('Coreg*TOF*MRA*'); %% change this
CoregInputs.RefImgPath=[pathOfNifti,filesep,fixedImg.name];
CoregInputs.Prefix='Coreg_';
Coregistration_job(CoregInputs);
disp('Coregistered TOF-MRA (fixed) and T1 MR (moving) co-registered!'); 
coregBrainMask=dir('Coreg_Brain_Mask_1.nii');
movefile(coregBrainMask.name,'Brain_Mask_1.nii');
coregTOFMRfile=dir('Coreg_TOF_MRA.nii');
movefile(coregTOFMRfile.name,'TOF_MRA_1.nii');

%% Perform coregistration between the MR navigators

% Nav 1+k (where k>0) is always the reference volume (stationary) and 
% Nav 1 is always the source (moving) volume. The transformation matrices
% obtained from this coregistration will be applied to the nifti mr Mask,
% to achieve spatial correspondence between the mr mask and the PET
% volumes.
clear CoregInputs
CoregInputs.Interp=0;
cd(pathOfNifti)
SourceMRnavFile=dir('*MRnav_1_*.nii')
CoregInputs.SourceImgPath=[pathOfNifti,filesep,SourceMRnavFile.name];
cd(pathOfNifti)
mrMaskFile=dir('MR*mask*1*.nii')
brainMaskFile=dir('Brain_Mask_1.nii');
TOFMRfile=dir('TOF_MRA_1.nii');
CoregInputs.MaskImgPath{1,:}=[pathOfNifti,filesep,mrMaskFile.name];
CoregInputs.MaskImgPath{2,:}=[pathOfNifti,filesep,brainMaskFile.name];
CoregInputs.MaskImgPath{3,:}=[pathOfNifti,filesep,TOFMRfile.name];
MRnavigators=dir('*MRnav*PET.nii');
for lp=1:length(MRnavigators)
    mrNavToSort{lp,:}=MRnavigators(lp).name;
end
sortedMRnav=natsort(mrNavToSort);
for lp=2:length(sortedMRnav)
    CoregInputs.Prefix=['Motion_',num2str(lp),'_'];
    CoregInputs.RefImgPath=[pathOfNifti,filesep,sortedMRnav{lp,:}]
    Coregistration_job(CoregInputs);
end

%% Move brain masks and petrous mask to seperate folder inside the moco
cd(pathOfNifti)
mkdir('MoCo_Brain_masks');
mkdir('MoCo_MR_masks');
mkdir('MoCo_TOF_MRA');
pathOfBrainMasksMRMoCo=[pathOfNifti,filesep,'MoCo_Brain_masks'];
pathOfMRmasksMRMoCo=[pathOfNifti,filesep,'MoCo_MR_masks'];
pathOfTOFMRAMRMoCo=[pathOfNifti,filesep,'MoCo_TOF_MRA'];
movefile('*Brain_Mask_*',pathOfBrainMasksMRMoCo);
movefile('*MR_mask_*',pathOfMRmasksMRMoCo);
movefile('*TOF_MRA_*',pathOfTOFMRAMRMoCo);
cd(pathOfBrainMasksMRMoCo)
motionFiles=dir('Motion*Brain*mask*');
for lp=1:length(motionFiles)
    brainFilesToSort{lp,:}=motionFiles(lp).name;
end
sortedBrains=natsort(brainFilesToSort);
for lp=1:length(sortedBrains)
    movefile(sortedBrains{lp,:},['Brain_mask_',num2str(lp+1),'.nii']);
end
cd(pathOfMRmasksMRMoCo)
motionFiles=dir('Motion*MR*mask*');
for lp=1:length(motionFiles)
    mrMaskFilesToSort{lp,:}=motionFiles(lp).name;
end
sortedMRmasks=natsort(mrMaskFilesToSort);
for lp=1:length(sortedMRmasks)
    movefile(sortedMRmasks{lp,:},['MR_mask_',num2str(lp+1),'.nii']);
end
cd(pathOfTOFMRAMRMoCo)
motionFiles=dir('Motion*TOF*MRA*');
for lp=1:length(motionFiles)
    tofMRAFilesToSort{lp,:}=motionFiles(lp).name;
end
sortedTOFMRA=natsort(tofMRAFilesToSort);
for lp=1:length(sortedTOFMRA)
    movefile(sortedTOFMRA{lp,:},['TOF_MRA_',num2str(lp+1),'.nii']);
end
cd(pathOfNifti)
mkdir('T1_MR');
pathOfT1mr=[pathOfNifti,filesep,'T1_MR'];
movefile('*T1_MR.*',pathOfT1mr);

cd(pathOfNifti)
mkdir('MR MoCo')
pathOfMRMoCoFiles=[pwd,filesep,'MR MoCo'];
movefile(pathOfBrainMasksMRMoCo,pathOfMRMoCoFiles)
movefile(pathOfMRmasksMRMoCo,pathOfMRMoCoFiles)
movefile(pathOfTOFMRAMRMoCo,pathOfMRMoCoFiles)
disp('MR based motion correction done...!');
pathOfBrainMasksMRMoCo=[pathOfMRMoCoFiles,filesep,'MoCo_Brain_masks'];
pathOfMRmasksMRMoCo=[pathOfMRMoCoFiles,filesep,'MoCo_MR_masks'];
pathOfTOFMRAMRMoCo=[pathOfMRMoCoFiles,filesep,'MoCo_TOF_MRA'];

%% Performing PET based MoCo only for the last ten frames. 

% Creation of folder structures
cd(pathOfNifti);
mkdir('PET MoCo');
pathOfPETMoCoFiles=[pwd,filesep,'PET MoCo'];
cd(pathOfPETMoCoFiles);
mkdir('MoCo_Brain_masks');
mkdir('MoCo_MR_masks');
mkdir('MoCo_TOF_MRA');
pathOfBrainMasksPETMoCo=[pathOfPETMoCoFiles,filesep,'MoCo_Brain_masks'];
pathOfMRmasksPETMoCo=[pathOfPETMoCoFiles,filesep,'MoCo_MR_masks'];
pathOfTOFMRAPETMoCo=[pathOfPETMoCoFiles,filesep,'MoCo_TOF_MRA'];
% Coregistration modules performing the PET motion correction.
clear CoregInputs;
CoregInputs.Interp=0; % nearest neighbor interpolation, as its a binary mask
CoregInputs.Prefix='Coreg_';
CoregInputs.SourceImgPath=[pathOfNifti,filesep,'TOF_MRA.nii'];
% cd(mrMoCoInputs.pathOfDICOMMRmask)
% convertDicomtoNii(cd,cd);
% oldFileName=dir('*.nii');
% movefile(oldFileName.name,'Native_MR_mask.nii');
% movefile('Native_MR_mask.nii',pathOfMRmasksMRMoCo);
CoregInputs.MaskImgPath{1,:}=[pathOfNifti,filesep,'MR_mask.nii'];
CoregInputs.MaskImgPath{2,:}=[pathOfNifti,filesep,'Brain_Mask.nii'];
cd(mrMoCoInputs.pathOfNiftiPET) % go inside the pet folder 
numberOfTotalPETfiles=length(dir('*.nii'));

for lp=(numberOfTotalPETfiles)-6:numberOfTotalPETfiles
    cd(mrMoCoInputs.pathOfNiftiPET);
    fileOfInterest=dir(['*',num2str(lp),'.nii']);
    CoregInputs.RefImgPath=[mrMoCoInputs.pathOfNiftiPET,filesep,fileOfInterest.name];
    Coregistration_job(CoregInputs);
    cd(pathOfNifti)
    oldFileName=dir('Coreg*Brain*');
    newFileName=['Brain_Mask_',num2str(lp),'.nii'];
    movefile(oldFileName.name,pathOfBrainMasksPETMoCo);
    cd(pathOfBrainMasksPETMoCo);
    movefile(oldFileName.name,newFileName);
    cd(pathOfNifti)
    oldFileName=dir('Coreg*MR*mask*');
    newFileName=['MR_mask_',num2str(lp),'.nii'];
    movefile(oldFileName.name,pathOfMRmasksPETMoCo);
    cd(pathOfMRmasksPETMoCo)
    movefile(oldFileName.name,newFileName);
    cd(pathOfNifti)
    oldFileName=dir('Coreg*TOF*MRA*');
    newFileName=['TOF_MRA_',num2str(lp),'.nii'];
    movefile(oldFileName.name,pathOfTOFMRAPETMoCo);
    cd(pathOfTOFMRAPETMoCo)
    movefile(oldFileName.name,newFileName);
end

disp('PET based MoCo performed for last 10 frames, as a back up measure...')
disp(['PET based MoCo files can be found in ',pathOfPETMoCoFiles,'!']);
    