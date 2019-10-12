%% This function is specifically designed for the IDIF protocol in AKH and will rotate the attenuation maps based on the transformation matrices obtained from the MR Navigators. 
%  Not for outside use - modifications "WILL" be needed.
%   
%  Author: Lalith Kumar Shiyam Sundar, M.Sc.
%  Date: Mar 23rd, 2018
%  Quantitative Imaging and Medical Physics, Medical University of Vienna.
%
%  Inputs: 
%       iRAM.pathOfStaticAC= Path where different static AC maps are present         
%       iRAM.pathOfMRnavigators = path where the MR navigators are present,
%  MR navigators should not be in mosaic format.
%
%  Outputs: 
%       The output is going to be a folder called "Dynamic AC maps"
%%------------------------------------------------------------------------%
%                               Program start
%%------------------------------------------------------------------------%

function []=idifRotateACmaps(iRAM)

%% Initialization and hardcoded variables 
cd(iRAM.pathOfStaticAC); 
cd ..
pathOfAC=cd;
outputFolder='Dynamic AC'; % Folder where the rotated attenuation maps will be stored.
workingFolder='Working area'; % buffer folder where all the temporary variables are stored.
mkdir(outputFolder); % create the output folder which contains the dynamic attenuation maps.
mkdir(workingFolder); % create the buffer folder.
pathToOutputFolder=[pathOfAC,filesep,outputFolder]; % paths to the newly created folder
pathToWorkingFolder=[pathOfAC,filesep,workingFolder]; % paths to the newly created folder 
cd ..
processedDataPath=cd;
%% Read the different type of attenuation maps

cd(iRAM.pathOfStaticAC);
ACnames=dir;
ACnames=ACnames(arrayfun(@(x) x.name(1), ACnames) ~= '.'); % read volunteer folders

%% Create the dynamic series folders in output folder.

cd(pathToOutputFolder); % go to the user defined path
for lp=1:length(ACnames)
    MoCoSeries{lp}=[ACnames(lp).name];
    mkdir(MoCoSeries{lp});
end

%% Create the working area for each attenuation map 

cd(pathToWorkingFolder); 
for lp=1:length(ACnames)
    workingArea{lp}=[ACnames(lp).name];
    mkdir(workingArea{lp});
end

%% Conversion of the low-dose CT to CT u-map using carney bilinear scaling.

cd(iRAM.pathOfStaticAC);
for lp=1:length(ACnames)
    if strcmp(ACnames(lp).name,'CT')
        lowDoseCTpath=[iRAM.pathOfStaticAC,filesep,ACnames(lp).name]; % find the low dose CT folder and save its path.
    else
    end
    if strcmp(ACnames(lp).name, 'DIXON')
        dixonSeriesPath=[iRAM.pathOfStaticAC,filesep,ACnames(lp).name]; % find the path of DIXON.
    else
    end
end

mkdir('CBS CT uMap'); % create a folder for the carney bilinear scaled CT umap;
CTuMapPath=[iRAM.pathOfStaticAC,filesep,'CBS CT uMap']; % path to the CT uMap which was scaled using carney bilinear scaling.

% setting up the hard-coded names for the generated image
lowDoseCTfile='Low_Dose_CT.nii';
DixonCompositeFile='Dixon-Composite-Image.nii';

% Read the dixon uMap and copy it to the reference folder of the CT series 
cd(dixonSeriesPath);
dixonUmapFolder=dir('*UMAP*');
dixonUmapPath=[dixonSeriesPath,filesep,dixonUmapFolder.name];
cd(dixonUmapPath);
dixonUmap=readImages(cd);
copyfile(dixonUmapPath,CTuMapPath); % folder where the Dixon dicom header information is stored
cd(dixonSeriesPath);
fNames=dir;
fNames=fNames(arrayfun(@(x) x.name(1), fNames) ~= '.'); % read volunteer folders
iterVar=0; %intialization
for lp=1:length(fNames)
    if isempty(strfind(fNames(lp).name,'UMAP'))
        iterVar=iterVar+1;
        DixonSeries{iterVar}=fNames(lp).name;
    else
    end
end

% convert the dixon series : Fat, in, opp, water to nifty files.

cd(pathToWorkingFolder);
for lp=1:length(workingArea)
    if strcmp(workingArea{lp},'DIXON')
        pathOfNiftiDixon=[pathToWorkingFolder,filesep,workingArea{lp}];
    end
    if strcmp(workingArea{lp},'CT')
        pathOfNiftiCT=[pathToWorkingFolder,filesep,workingArea{lp}];
    end
end
for lp=1:length(DixonSeries)
    tempDixonPath=[dixonSeriesPath,filesep,DixonSeries{lp}];
    cd(tempDixonPath);
    convertDicomtoNii(cd,cd);
    niftyFile=dir('*.nii');
    movefile(niftyFile.name,[DixonSeries{lp},'.nii']);
    movefile([DixonSeries{lp},'.nii'],pathOfNiftiDixon);
end

% Preparing a Dixon composite image using the fat, water, inphase, outphase maps
 
cd(pathOfNiftiDixon); % go to the folder containing the individual dixon volumes (nifti format)
NiftiDcmFiles=dir;NiftiDcmFiles=NiftiDcmFiles(arrayfun(@(x) x.name(1), NiftiDcmFiles) ~= '.'); % reading the files inside the folder.
DixonComposite=zeros(size(spm_read_vols(spm_vol(NiftiDcmFiles(1).name)))); % Initialize an empty volume grid.
for lp=1:length(NiftiDcmFiles)
    DixonComposite=DixonComposite+spm_read_vols(spm_vol(NiftiDcmFiles(lp).name));
end
NiftiDixonHdr=spm_vol(NiftiDcmFiles(1).name);
NiftiDixonHdr.fname='Dixon-Composite-Image.nii';
spm_write_vol(NiftiDixonHdr,DixonComposite);

%% Preparing the low-dose CT image 

cd(lowDoseCTpath);
convertDicomtoNii(cd,cd);
niftiFile=dir('*.nii');
movefile(niftiFile.name,lowDoseCTfile);
movefile(lowDoseCTfile,pathOfNiftiCT);
cd(pathOfNiftiCT);

% Segmenting the bed from the CT image 

OrgCT=spm_read_vols(spm_vol(lowDoseCTfile)); % reads in the Nifti CT images.
figure,imshow3D(OrgCT ); % Displaying the mid axial slice of the CT volume
title('Original CT image with bed'); 

% 2.) Performing a otsu thresholding.

levelThresh=multithresh(OrgCT,1); % segmenting the foreground and background
SkullThresh=imquantize(OrgCT,levelThresh); %labelling the objects.
foreGround=(SkullThresh==2); % segmenting the bright object - foreground.
[BinaryCT,~]=VolumePartitioner(foreGround); % getting only the head.

for lp=1:size(BinaryCT,3)
    OnlyTheHead(:,:,lp)=imfill(BinaryCT(:,:,lp),'holes'); % filling up the holes in the head, this is done to segment out only the head
end
OnlyTheHead=OnlyTheHead>0; % converting into logical;
figure,imshow3D(OnlyTheHead);title('Binary mask with the holes in the head filled');
tempImg=OnlyTheHead.*OrgCT; 

% Scale the CT head, so that the hounsfield units are in the range of 0 to 4000
scaledImg=tempImg+1000;
OnlyTheHead=scaledImg.*OnlyTheHead; % using the head mask to get only the head from the scaled Image.

% Stealing the NIFTI header and writing the CT volume without the bed as NIFTI file
% We need a NIFTI header for the CT volume without the bed, which we have 
% generated. So we are gonna use the NIFTI header from the original CT volume.

NiftiHdrCT=spm_vol(lowDoseCTfile); % read the nifti header of the original CT file using the 'spm_vol' function.
NiftiHdrCT.fname='CTwithoutBed.nii'; % name for the CT volume without the bed - i have hardcoded it, you can change it.
spm_write_vol(NiftiHdrCT,OnlyTheHead); % writing the CT volume without the bed as a NIFTI file.

% Code-patch for setting the origin to the center of the image volume, instead of
% manual clicks.

st.vol=spm_vol('CTwithoutBed.nii');
vs=st.vol.mat\eye(4);
vs(1:3,4)=(st.vol.dim+1)/2;
spm_get_space(st.vol.fname,inv(vs));


%% Perform co-registration between the dixon and CT without the bed -> Transforming from PET/CT to PET/MR

CoregInputs.SourceImgPath=[pathOfNiftiCT,filesep,'CTwithoutBed.nii'];
CoregInputs.RefImgPath=[pathOfNiftiDixon,filesep,DixonCompositeFile];
CoregInputs.MaskImgPath={''};
CoregInputs.Prefix='Coreg_';
CoregInputs.Interp=1; % was 4th degree spline '4' before, now its trilinear interpolation. 
Coregistration_job(CoregInputs);

%% Cleaning the remaining bed, if there is any.

CTimgToScale=spm_read_vols(spm_vol('Coreg_CTwithoutBed.nii'));


%% Perform bilinear scaling and push the CT-Umap to the dixon header (reference CT u-Map)

CTuMap=carneyBilinearScaling(lowDoseCTpath,CTimgToScale);

%% Making sure that CT-umap doesn't have the bed.

[CTuMapOnly,~]=VolumePartitioner(CTuMap);
for lp=1:size(CTuMapOnly,3)
    headMask(:,:,lp)=imfill(CTuMapOnly(:,:,lp),'holes'); % filling up the holes in the head, this is done to segment out only the head
end
headMask=headMask>0; % converting into logical;
figure,imshow3D(headMask);title('Binary mask with the holes in the head filled');
CTuMapCleaned=CTuMap.*headMask; 
%pause; % wait for keypress, a keypress indicates, that the u-map is free of bed and the program continues.

pushDataToDicom(CTuMapPath,CTuMapCleaned);

% move the low dose CT from AC maps folder to the processed data folder 
%pause(5)
movefile (lowDoseCTpath,processedDataPath)
cd(iRAM.pathOfStaticAC);
movefile('CBS CT uMap','CT');


%% Handle RESOLUTE maps, resolute maps need to be in alignment with the Dixon composite image.
% pathToNiftiResolute=[pathToWorkingFolder,filesep,'RESOLUTE'];
% cd([iRAM.pathOfStaticAC,filesep,'RESOLUTE'])
% convertDicomtoNii(cd,cd);
% oldResoluteFile=dir('*nii');
% movefile(oldResoluteFile.name,'RESOLUTE.nii')
% movefile('RESOLUTE.nii',pathToNiftiResolute)
% cd(pathToNiftiResolute)
% resoluteFile=dir('*nii');
% clear CoregInputs;
% CoregInputs.SourceImgPath=[pathToNiftiResolute,filesep,resoluteFile.name];
% CoregInputs.RefImgPath=[pathOfNiftiDixon,filesep,DixonCompositeFile];
% CoregInputs.Interp=1;
% CoregInputs.MaskImgPath={''};
% CoregInputs.Prefix='Coreg_'
% Coregistration_job(CoregInputs);
% delete 'RESOLUTE.nii';
% coreg_ResoluteFile=dir('Coreg*');
% movefile(coreg_ResoluteFile.name,'RESOLUTE.nii');
% resoluteVolume=spm_read_vols(spm_vol('RESOLUTE.nii'));
% pathOfDicomResolute=[iRAM.pathOfStaticAC,filesep,'RESOLUTE'];
% cd(pathOfDicomResolute)
% delete *.IMA
% copyfile(dixonUmapPath,pathOfDicomResolute)
% pushDataToDicom(pathOfDicomResolute,resoluteVolume); 

%% 
% cd(iRAM.pathOfStaticAC)
% rmdir 'DIXON' s  
% cd(pathToOutputFolder)
% rmdir 'DIXON' s
% cd(pathToWorkingFolder)
% rmdir 'DIXON' s

%% Find the number of MR navigators and generate equal number of CT u-maps 

cd(iRAM.pathOfMRnavigators)
subFolders=dir;
subFolders=subFolders(arrayfun(@(x) x.name(1), subFolders) ~= '.'); % read volunteer folders
numberOfMRnavigators=length(subFolders);
cd(iRAM.pathOfStaticAC);
ACnames=dir;
ACnames=ACnames(arrayfun(@(x) x.name(1), ACnames) ~= '.');

% Copy the static folders to dynamic folders based on the number of mr  
% navigators present
cd(pathToOutputFolder);
wb=timebar('Preparing the buffer u-maps...');
for olp=1:length(ACnames)
    for ilp=1:(numberOfMRnavigators)
        if strcmp(ACnames(olp).name,'DIXON')
            staticPath=dixonUmapPath; 
            continue
        else
            staticPath=[iRAM.pathOfStaticAC,filesep,ACnames(olp).name];
        end
        ACseriesName=[ACnames(olp).name,'_uMap_',num2str(ilp)];
        folderOfInterest= dir([ACnames(olp).name,'*']);
        destinationPath=[pathToOutputFolder,filesep,folderOfInterest.name,filesep,ACseriesName];
        copyfile (staticPath,destinationPath)
    end
    timebar(wb,olp/length(ACnames));
end
close(wb)

%% Convert the static AC maps to nifti and store it in working folder for further processing.

cd(iRAM.pathOfStaticAC);
for lp=1:length(ACnames)
    if strcmp(ACnames(lp).name,'DIXON')
        tempPath=dixonUmapPath;
    else
        tempPath=([iRAM.pathOfStaticAC,filesep,ACnames(lp).name]);
    end
    convertDicomtoNii(tempPath,tempPath);
    whereToMove=[pathToWorkingFolder,filesep,ACnames(lp).name];
    niftyFile=dir('*nii');
    movefile (niftyFile.name,[ACnames(lp).name,'.nii']);
    movefile([ACnames(lp).name,'.nii'],whereToMove);
end


%% Transform the MR navigators to world coordinates, also transform the MR navigators to CT-space.

cd(lowDoseCTpath);
CT=readImages(cd);
cd(iRAM.pathOfMRnavigators);
subFolders=dir;
subFolders=subFolders(arrayfun(@(x) x.name(1), subFolders) ~= '.'); % read volunteer folders

% sort the mr navigator folders using "natsort" 

for lp=1:length(subFolders)
    foldersToSort{lp,:}=subFolders(lp).name;
end

sortedFolders=natsort(foldersToSort);
wb=timebar('Transforming images from voxel space to world-space (mm)...');

for lp=1:length(foldersToSort)
    cd([iRAM.pathOfMRnavigators,filesep,sortedFolders{lp}])
    MRnavigator=readImages(cd,'align');
    outMRnav=transformVoxelSpace(MRnavigator,CT);
    cd(iRAM.pathOfMRnavigators);
    writeDicomImage(outMRnav,['Transform_MRnav_Dixon_',num2str(lp)]);
    timebar(wb,lp/length(foldersToSort));
end
close(wb);
disp('Transforming MR navigators to DIXON space...');
%% Convert the DICOM MR navigators to Nifti MR navigators.

cd(iRAM.pathOfMRnavigators);
mkdir('Nifty_MRnav_Dixon');
pathOfNiftyMRnav=[iRAM.pathOfMRnavigators,filesep,'Nifty_MRnav_Dixon'];


% Convert the transformed MR navigators to nifti files using SPM

transMRnavFiles=dir('Transform*Dixon*');
for lp=1:length(transMRnavFiles)
    pathOfDicom=[iRAM.pathOfMRnavigators,filesep,transMRnavFiles(lp).name];
    convertDicomtoNii(pathOfDicom,pathOfDicom);
    cd(pathOfDicom)
    niftyFile=dir('*.nii');
    movefile(niftyFile.name,[transMRnavFiles(lp).name,'.nii'])
    movefile([transMRnavFiles(lp).name,'.nii'],pathOfNiftyMRnav)
end


%% Perform rigid registration using the MR navigators and apply the transformation matrices to the umaps

% Since RESOLUTE and UTE are going to have different matrix sizes, they are
% handled separately. Here other u-maps are handled.

iter=0;
for lp=1:length(ACnames)
    tempPath=[pathToWorkingFolder,filesep,ACnames(lp).name,filesep,[ACnames(lp).name,'.nii']];
    iter=iter+1;
    otherImages{iter,:}=tempPath;
    filesInSameSpace{iter,:}=ACnames(lp).name;   
end
CoregInputs.MaskImgPath=otherImages;


%%
CoregInputs.Interp=1; % it was 4th degree spline interpolation ('4'-spm trigger) before. now it's trilinear interpolation '1'
cd(pathOfNiftyMRnav);
SourceMRnavFile=dir('*_1.nii'); % the image thats gonna be jiggled.
CoregInputs.SourceImgPath=[pathOfNiftyMRnav,filesep,SourceMRnavFile.name];
cd(pathOfNiftyMRnav)
for lp=2:numberOfMRnavigators
    CoregInputs.Prefix=['First_To_',num2str(lp),'_'];
    RefMRnavFile=dir(['*_',num2str(lp),'.nii']); % the image which is gonna stay put.
    CoregInputs.RefImgPath=[pathOfNiftyMRnav,filesep,RefMRnavFile.name];
    Coregistration_job(CoregInputs)
end



%% Push the transformed u-map data to the respective DICOM counter parts.

  wb=timebar('Pushing rotated binary data to their respective DICOM AC maps...');
for olp=1:length(filesInSameSpace)
    for ilp=2:numberOfMRnavigators
        tempPath=[pathToWorkingFolder,filesep,filesInSameSpace{olp}];
        cd(tempPath);
        fileName=(['*','_',num2str(ilp),'_*']);
        fileToRead=dir(fileName);
        disp(['Pushing ',fileToRead.name,' now!']);
        volumeToPush=spm_read_vols(spm_vol(fileToRead.name));
        cd([pathToOutputFolder,filesep,filesInSameSpace{olp}]);
        folderToPush=dir(['*_uMap_',num2str(ilp)]);
        pathToPush=[pathToOutputFolder,filesep,filesInSameSpace{olp},filesep,folderToPush.name];
        pushDataToDicom(pathToPush,volumeToPush);
    end
timebar(wb,olp/length(ACnames));

end
close(wb)

%% Handling the UTE space attenuation maps 

% coregister the RESOLUTE attenuation map to 

% 
% % Transform the MR navigators to world coordinates, also transform the MR navigators to CT-space.
% pathOfUTE=[iRAM.pathOfStaticAC,filesep,'UTE'];
% cd(pathOfUTE);
% UTE=readImages(cd);
% sortedFolders=natsort(foldersToSort);
% wb=timebar('Transforming images from voxel space to world-space (mm)...');
% 
% for lp=1:length(foldersToSort)
%     cd([iRAM.pathOfMRnavigators,filesep,sortedFolders{lp}])
%     MRnavigator=readImages(cd,'align');
%     outMRnav=transformVoxelSpace(MRnavigator,UTE);
%     cd(iRAM.pathOfMRnavigators);
%     writeDicomImage(outMRnav,['Transform_MRnav_UTE_',num2str(lp)]);
%     timebar(wb,lp/length(foldersToSort));
% end
% close(wb);
% disp('Transforming MR navigators to UTE space...');
% 
% %% Convert the DICOM MR navigators to Nifti MR navigators.
% 
% cd(iRAM.pathOfMRnavigators);
% mkdir('Nifty_MRnav_UTE');
% pathOfNiftyMRnav=[iRAM.pathOfMRnavigators,filesep,'Nifty_MRnav_UTE'];
% 
% 
% % Convert the transformed MR navigators to nifti files using SPM
% 
% transMRnavFiles=dir('Transform*UTE*');
% for lp=1:length(transMRnavFiles)
%     pathOfDicom=[iRAM.pathOfMRnavigators,filesep,transMRnavFiles(lp).name];
%     convertDicomtoNii(pathOfDicom,pathOfDicom);
%     cd(pathOfDicom)
%     niftyFile=dir('*.nii');
%     movefile(niftyFile.name,[transMRnavFiles(lp).name,'.nii'])
%     movefile([transMRnavFiles(lp).name,'.nii'],pathOfNiftyMRnav)
% end
% 
% 
% %% Perform rigid registration using the MR navigators and apply the transformation matrices to the umaps
% 
% % Since RESOLUTE and UTE are going to have different matrix sizes, they are
% % handled here.
% 
% clear tempPath filesInSameSpace otherImages CoregInputs
% iter=0;
% for lp=1:length(ACnames)
%     if strcmp(ACnames(lp).name,'RESOLUTE') || strcmp(ACnames(lp).name,'UTE')
%         tempPath=[pathToWorkingFolder,filesep,ACnames(lp).name,filesep,[ACnames(lp).name,'.nii']];
%         iter=iter+1;
%         otherImages{iter,:}=tempPath;
%         filesInSameSpace{iter,:}=ACnames(lp).name;
%     else
%         continue
%     end
% end
% CoregInputs.MaskImgPath=otherImages;
% 
% %% running the coregistration
% 
% CoregInputs.Interp=1; % trilinear interpolation 
% cd(pathOfNiftyMRnav);
% SourceMRnavFile=dir('*_1.nii'); % the image thats gonna be jiggled.
% CoregInputs.SourceImgPath=[pathOfNiftyMRnav,filesep,SourceMRnavFile.name];
% cd(pathOfNiftyMRnav)
% for lp=2:numberOfMRnavigators
%     CoregInputs.Prefix=['First_To_',num2str(lp),'_'];
%     RefMRnavFile=dir(['*_',num2str(lp),'.nii']); % the image which is gonna stay put.
%     CoregInputs.RefImgPath=[pathOfNiftyMRnav,filesep,RefMRnavFile.name];
%     Coregistration_job(CoregInputs)
% end
% %% Push the transformed u-map data to the respective DICOM counter parts.
% 
% wb=timebar('Pushing rotated binary data to their respective DICOM AC maps...');
% for olp=1:length(filesInSameSpace)
%     if strcmp(filesInSameSpace{olp,:},'RESOLUTE')
%         startPoint=1;
%         cd([pathToWorkingFolder,filesep,filesInSameSpace{olp}])
%         oldFile=dir('RESOLUTE.nii');
%         movefile(oldFile.name,'First_To_1_RESOLUTE.nii');
%     else
%         if strcmp(filesInSameSpace{olp,:},'UTE')
%             startPoint=2;
%         end
%     end
%     for ilp=startPoint:numberOfMRnavigators
%         tempPath=[pathToWorkingFolder,filesep,filesInSameSpace{olp}];
%         cd(tempPath);
%         fileName=(['*','_',num2str(ilp),'_*']);
%         fileToRead=dir(fileName);
%         disp(['Pushing ',fileToRead.name,' now!']);
%         volumeToPush=spm_read_vols(spm_vol(fileToRead.name));
%         cd([pathToOutputFolder,filesep,filesInSameSpace{olp}]);
%         folderToPush=dir(['*_uMap_',num2str(ilp)]);
%         pathToPush=[pathToOutputFolder,filesep,filesInSameSpace{olp},filesep,folderToPush.name];
%         pushDataToDicom(pathToPush,volumeToPush);
%     end
% timebar(wb,olp/length(ACnames));
% 
% end
% close(wb)
% 
% %% Creating a histogram overlay for visualizing the attenuation distribution for each map. 
% 
% cd(pathToOutputFolder)% go to the dynamic folder
% acFolders=dir;
% acFolders=acFolders(arrayfun(@(x) x.name(1), acFolders) ~= '.'); % read volunteer folders
% for lp=1:length(acFolders)
%     temp=[pathToOutputFolder,filesep,acFolders(lp).name];
%     cd(temp);
%     subFolders=dir;
%     subFolders=subFolders(arrayfun(@(x) x.name(1), subFolders) ~= '.'); % read volunteer folders
%     figure,
%     for ilp=1:length(subFolders)
%         cd([temp,filesep,subFolders(ilp).name]);
%         imgVol=MPRage(cd);
%         histogram(imgVol(:));
%         hold on, title(subFolders(ilp).name);
%         cd .. 
%     end
% end
    
    