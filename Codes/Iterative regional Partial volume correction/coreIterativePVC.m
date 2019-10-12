function [pvcOutputs] = coreIterativePVC(PVC,r)

%% Hard-coded variables : 
numberOfThresholds=10;                                                      % Define the number of thresholds for otsu thresholding (the number 20 is basically chosen in a random manner)
widthOfSE=10;
cd(PVC.whereToProcess)
mkdir('PVC')
workingFolder=[PVC.whereToProcess,filesep,'PVC'];
cd(workingFolder);
mkdir('PVC masks');
pathToPVCMasks=[workingFolder,filesep,'PVC masks'];
cd(pathToPVCMasks);
mkdir('Artery masks');
pathToArteryMasks=[pathToPVCMasks,filesep,'Artery masks'];
mkdir('Brain masks');
pathToBrainMasks=[pathToPVCMasks,filesep,'Brain masks'];

%% Read in PET volumes 

cd(PVC.pathOfPET); % go to the folder containing the nifti pet volumes 
PETfiles=dir; PETfiles=PETfiles(arrayfun(@(x) x.name(1), PETfiles) ~= '.'); % read the PET files in the folder
parfor lp=1:length(PETfiles)
    PETfilesToSort{lp,:}=PETfiles(lp).name;
end
sortedPETfiles=natsort(PETfilesToSort); % sorting the pet files based on their name.
parfor lp=1:length(PETfiles)
    PET{lp}=spm_read_vols(spm_vol(sortedPETfiles{lp}))./(PVC.crossCalibrationFactor); % loading the pet files.
    petHdr{lp}=spm_vol(sortedPETfiles{lp});
    disp(['Reading ',sortedPETfiles{lp},'...']);
end
PETframeStruct=petHdr{1};

%% Normalized psfSigma based on PET voxel size. 
xDim=abs(PETframeStruct.mat(1,1)); % get the voxel size (x dimension)
yDim=abs(PETframeStruct.mat(2,2));% get the voxel size (y dimension)
zDim=abs(PETframeStruct.mat(3,3));% get the voxel size (z dimension)
voxelSize=[xDim yDim zDim];

% normalize the FWHM to voxel size 

if numel(PVC.psfFWHM)==1
    psfFWHM=[PVC.psfFWHM PVC.psfFWHM PVC.psfFWHM]; % isotropic filtering.                
else
end

psfFWHM  = psfFWHM./voxelSize; % normalizing the fwhm to the pet voxel dimensions.
psfSigma = psfFWHM/sqrt(8*log(2)); % calculating the sigma from fwhm    

%% calculate pet mid time

[PETframingMidPoint]=getPETmidTime();

%% Read in MR masks based on MR based MoCo.

cd(PVC.pathToMRcorrespondenceFile)
mrCorrespondenceFile=dir('MR_correspondence.txt');
mrCorrespondence=readtable(mrCorrespondenceFile.name);
mrCorrespondence=mrCorrespondence.Var1;

%% Generate a folder which contains MR masks for each pet frame : correspondence achieved using "MR correspondence file"

parfor lp=1:length(mrCorrespondence)
    cd(PVC.pathOfhighResMask)
    fileToRead=dir(['MR_mask_',num2str(mrCorrespondence(lp)),'.nii']); 
    volToWrite=spm_read_vols(spm_vol(fileToRead.name));
    dummyHdr=(spm_vol(fileToRead.name));
    dummyHdr.fname=['Artery_mask_',num2str(lp),'.nii'];
    cd(pathToArteryMasks);
    spm_write_vol(dummyHdr,volToWrite);
    disp(['Writing ',dummyHdr.fname,'!'])
    cd(PVC.pathOfBrainMasks);
    fileToRead=dir(['Brain_mask_',num2str(mrCorrespondence(lp)),'.nii']);
    volToWrite=spm_read_vols(spm_vol(fileToRead.name));
    dummyHdr=spm_vol(fileToRead.name);
    dummyHdr.fname=['Brain_mask_',num2str(lp),'.nii'];
    cd(pathToBrainMasks);
    spm_write_vol(dummyHdr,volToWrite);
    disp(['Writing ',dummyHdr.fname,'!']);
    
end

 
% cd(PVC.pathOfpMoCoBrainMasks)
% mrFiles=dir; mrFiles=mrFiles(arrayfun(@(x) x.name(1), mrFiles) ~= '.'); % read the PET files in the folder
% for lp=1:length(mrFiles)    
%     cd(PVC.pathOfpMoCoBrainMasks)
%     Num=regexp(mrFiles(lp).name,'\d');
%     indexToExtract(lp)=str2num(mrFiles(lp).name(Num));
%     FOI=dir(['*_',num2str(indexToExtract(lp)),'.nii']);
%     if isempty(FOI)
%        folderName=['B_',num2str(indexToExtract(lp))];
%        fileToMove=[pwd,filesep,folderName,'*'];
%     else
%        fileToMove=[pwd,filesep,FOI.name];
%     end
%     cd(pathToBrainMasks)
%     fileToDelete=dir(['*',num2str(indexToExtract(lp)),'.nii']);
%     delete(fileToDelete.name);
%     copyfile(fileToMove,pathToBrainMasks);   
% end
% cd(PVC.pathOfpMoCohighResMask)
% mrFiles=dir; mrFiles=mrFiles(arrayfun(@(x) x.name(1), mrFiles) ~= '.'); % read the PET files in the folder
% for lp=1:length(mrFiles)    
%     cd(PVC.pathOfpMoCohighResMask)
%     Num=regexp(mrFiles(lp).name,'\d');
%     indexToExtract(lp)=str2num(mrFiles(lp).name(Num));
%     FOI=dir(['*_',num2str(indexToExtract(lp)),'.nii']);
%     if isempty(FOI)
%        folderName=['A_',num2str(indexToExtract(lp))];
%        fileToMove=[pwd,filesep,folderName,'*'];
%     else
%        fileToMove=[pwd,filesep,FOI.name];
%     end
%     cd(pathToArteryMasks)
%     fileToDelete=dir(['*',num2str(indexToExtract(lp)),'.nii']);
%     delete(fileToDelete.name);
%     copyfile(fileToMove,pathToArteryMasks); 
%     if isempty(FOI)
%     else
%         movefile(FOI.name,['Artery_mask_',num2str(indexToExtract(lp)),'.nii']);
%     end
% end

%% Read the artery and the brain masks in a single go.

for lp=1:length(mrCorrespondence);
    cd(pathToArteryMasks)
    FOI=dir(['*_',num2str(lp),'.nii']);
    if isempty(FOI)
        manuallyRotatedMRmask=dir(['A_',num2str(lp),'*']);
        cd([pathToArteryMasks,filesep,manuallyRotatedMRmask.name,filesep,'DCM000']);
        fileNames=dir('*.IMA');
        if isempty(fileNames)
            amideToMatlab(cd);
        end
        mrMaskInDicomOrientation=readImages(cd,'align');
        mrMaskInNiftiOrientation=rot90(mrMaskInDicomOrientation.volumes,-2);
        mrMaskInNiftiOrientation=double(mrMaskInNiftiOrientation>0);
        targetMRmask{lp}=mrMaskInNiftiOrientation;
        disp(['Reading ',manuallyRotatedMRmask.name,'...']);
    else
        targetMRmask{lp}=spm_read_vols(spm_vol(['Artery_mask_',num2str((lp)),'.nii']));
        if PVC.maskSwitch=='c' && lp < 26   % cervical segment
           [startInfo,endInfo]=StartEndImgInfo(targetMRmask{lp});
           maskOfInterest=((targetMRmask{lp}));
           %maskOfInterest(:,:,(endInfo-7):endInfo)=targetMRmask{lp}(:,:,(endInfo-7):endInfo);
           targetMRmask{lp}=targetMRmask{lp}.*maskOfInterest;
        end
        disp(['Reading ','Artery_mask_',num2str((lp)),'.nii']);
    end
    recoveryCoefficients{lp}=imgaussfilt3(double(targetMRmask{lp}),psfSigma);
    cd(pathToBrainMasks)
    FOI=dir(['*_',num2str((lp)),'.nii']);
    if isempty(FOI)
       manuallyRotatedBrainMask=dir(['B_',num2str(lp),'*']);
       cd([pathToBrainMasks,filesep,manuallyRotatedBrainMask.name,filesep,'DCM000'])
       fileNames=dir('*.IMA');
       if isempty(fileNames)
           amideToMatlab(cd);
       end
       brainMaskInDicomOrientation=readImages(cd,'align');
       brainMaskInNiftiOrientation=rot90(brainMaskInDicomOrientation.volumes,-2);
       brainMaskInNiftiOrientation=double(brainMaskInNiftiOrientation>0);
       brainMask{lp}=brainMaskInNiftiOrientation.*(~targetMRmask{lp});
       disp(['Reading ',manuallyRotatedBrainMask.name,'...']);
    else
       brainMask{lp}=spm_read_vols(spm_vol(['Brain_mask_',num2str(lp),'.nii']));
       disp(['Reading','Brain_mask_',num2str(lp),'.nii']);
   end
end

% % 
% for lp=1:length(mrCorrespondence)
%     cd(PVC.pathOfhighResMask)
%     fileToRead=dir(['MR_mask_',num2str(mrCorrespondence(lp)),'.nii']); 
%     dummyHdr=(spm_vol(fileToRead.name));
%     targetMRmask{lp}=(spm_read_vols(dummyHdr))>0;
%     manuallyRotatedMRmask=dir(['A_',num2str(mrCorrespondence(lp)),'*']);
%     if isempty(manuallyRotatedMRmask)
%     else
%         cd([PVC.pathOfhighResMask,filesep,manuallyRotatedMRmask.name,filesep,'DCM000']);
%         mrMaskInDicomOrientation=readImages(cd,'align');
%         mrMaskInNiftiOrientation=rot90(mrMaskInDicomOrientation.volumes,-2);
%         mrMaskInNiftiOrientation=double(mrMaskInNiftiOrientation>0);
%         targetMRmask{lp}=mrMaskInNiftiOrientation;
%         disp(['Reading ',manuallyRotatedMRmask.name,'...']);
%     end    
%     % Creation of recovery coefficients, the header info remains unchanged
%     % to preserve the values between 0 and 1
%     recoveryCoefficients{lp}=imgaussfilt3(double(targetMRmask{lp}),psfSigma);
%     cd(PVC.pathOfBrainMasks);
%     brainMaskFile=dir(['Brain_mask_',num2str(mrCorrespondence(lp)),'.nii']);
%     brainMask{lp}=spm_read_vols(spm_vol(brainMaskFile.name)).*(~targetMRmask{lp});
%     manuallyRotatedBrainMask=dir(['B_',num2str(mrCorrespondence(lp)),'*']);
%     if isempty(manuallyRotatedBrainMask)
%     else
%         cd([PVC.pathOfBrainMasks,filesep,manuallyRotatedBrainMask.name,filesep,'DCM000'])
%         brainMaskInDicomOrientation=readImages(cd,'align');
%         brainMaskInNiftiOrientation=rot90(brainMaskInDicomOrientation.volumes,-2);
%         brainMaskInNiftiOrientation=double(brainMaskInNiftiOrientation>0);
%         brainMask{lp}=brainMaskInNiftiOrientation.*(~targetMRmask{lp});
%         disp(['Reading ',manuallyRotatedBrainMask.name,'...']);
%     end     
%     disp(['Reading ',brainMaskFile.name,'...']);
% end
disp ('PET files read!');
disp ('High resolution MR mask files (brain+artery) read!');
disp ('Recovery coefficients calculated...!');

% %% Read in MR masks based on PET based MoCo, if pet based moco switch is set to true
% 
% if PVC.petMoCo==true
%    cd(PVC.pathOfpMoCohighResMask);
%    mrFiles=dir; mrFiles=mrFiles(arrayfun(@(x) x.name(1), mrFiles) ~= '.'); % read the PET files in the folder
% 
%    for lp=1:length(mrFiles)
%        cd(PVC.pathOfpMoCohighResMask);
%        Num=regexp(mrFiles(lp).name,'\d');
%        indexToExtract(lp)=str2num(mrFiles(lp).name(Num));
%        FOI=dir(['*',num2str(indexToExtract(lp)),'.nii']);
%        if isempty(FOI)
%            manuallyRotatedMRmask=dir(['A_',num2str(indexToExtract(lp)),'*']);
%            cd([PVC.pathOfpMoCohighResMask,filesep,manuallyRotatedMRmask.name,filesep,'DCM000']);
%            mrMaskInDicomOrientation=readImages(cd,'align');
%            mrMaskInNiftiOrientation=rot90(mrMaskInDicomOrientation.volumes,-2);
%            mrMaskInNiftiOrientation=double(mrMaskInNiftiOrientation>0);
%            targetMRmask{indexToExtract(lp)}=mrMaskInNiftiOrientation;
%            disp(['Reading ',manuallyRotatedMRmask.name,'...']);
%        else
%            targetMRmask{indexToExtract(lp)}=spm_read_vols(spm_vol(['MR_mask_',num2str(indexToExtract(lp)),'.nii']));
%            disp(['Reading ','MR_mask_',num2str(indexToExtract(lp)),'.nii']);
%        end
%        cd(PVC.pathOfpMoCoBrainMasks);
%        Num=regexp(mrFiles(lp).name,'\d');
%        indexToExtract(lp)=str2num(mrFiles(lp).name(Num));
%        FOI=dir(['*',num2str(indexToExtract(lp)),'.nii']);
%        if isempty(FOI)
%            manuallyRotatedBrainMask=dir(['B_',num2str(indexToExtract(lp)),'*']);
%            cd([PVC.pathOfpMoCoBrainMasks,filesep,manuallyRotatedBrainMask.name,filesep,'DCM000'])
%            brainMaskInDicomOrientation=readImages(cd,'align');
%            brainMaskInNiftiOrientation=rot90(brainMaskInDicomOrientation.volumes,-2);
%            brainMaskInNiftiOrientation=double(brainMaskInNiftiOrientation>0);
%            brainMask{indexToExtract(lp)}=brainMaskInNiftiOrientation.*(~targetMRmask{indexToExtract(lp)});
%            disp(['Reading ',manuallyRotatedBrainMask.name,'...']);
%        else
%            brainMask{indexToExtract(lp)}=spm_read_vols(spm_vol(['Brain_mask_',num2str(indexToExtract(lp)),'.nii']));
%            disp(['Reading','Brain_mask_',num2str(indexToExtract(lp)),'.nii']);
%        end
%       
%        recoveryCoefficients{indexToExtract(lp)}=imgaussfilt3(double(targetMRmask{indexToExtract(lp)}),psfSigma);
%        disp(['New recovery coefficiencts calculated for frames ',num2str(indexToExtract(lp)),'...']); 
%        
%    end
% else
% end


%% Create background mantel from the edges of the carotid artery. 

SE = strel('cube',widthOfSE);   
parfor lp=1:length(PETfiles)
    backgroundMantel{lp}=imdilate(targetMRmask{lp},SE)-targetMRmask{lp};
    disp(['Background mantel ',num2str(lp),' created!']);
    significantBrain{lp}=backgroundMantel{lp}.*brainMask{lp};
    %mixedZone{lp}=backgroundMantel{lp}-significantBrain{lp};
    brainSpillOutGeometry{lp}=(imgaussfilt3(double(significantBrain{lp}),psfSigma))>0.03; 
    arterySpillOutGeometry{lp}=recoveryCoefficients{lp}>0.03;
    cranialMixedZone{lp}=(brainSpillOutGeometry{lp}.*~significantBrain{lp});
    cranialMixedZone{lp}=cranialMixedZone{lp}.*(arterySpillOutGeometry{lp}.*~targetMRmask{lp});
    caudalMixedZone{lp}=arterySpillOutGeometry{lp}-(targetMRmask{lp}+cranialMixedZone{lp});
end
disp('Cranial Mixed zone geoemetry calculated...');
disp('Part of the brain which contributes to the spill-over to carotid calculated...');
disp('Background mantel calculated...');



%% Visualize the intensities in the above mentioned segments.
figure,
ylabel('Activity kBq/mL');xlabel('Time in seconds');title('Raw Activity analysis');
arteryProfile=animatedline('Color','red','Marker','o','LineWidth',3);
brainProfile=animatedline('Color','green','Marker','o','LineWidth',3);
cranialMixedZoneProfile=animatedline('Color','blue','Marker','o','LineWidth',3);
caudalMixedZoneProfile=animatedline('Color','black','Marker','o','LineWidth',3);
arteryProfile.DisplayName='Artery activity';
brainProfile.DisplayName='Brain activity';
cranialMixedZoneProfile.DisplayName='Cranial mixed zone activity';
caudalMixedZoneProfile.DisplayName='Caudal mixed zone activity';
for lp=1:length(significantBrain)
    brainTAC(lp)=median(PET{lp}(significantBrain{lp}>0))/1000;
    arteryTAC(lp)=mean(PET{lp}(targetMRmask{lp}>0))/1000;
    cranialTAC(lp)=median(PET{lp}(cranialMixedZone{lp}>0))/1000;
    caudalTAC(lp)=median(PET{lp}(caudalMixedZone{lp}>0))/1000;
    addpoints(brainProfile,PETframingMidPoint(lp),brainTAC(lp));
    hold on
    if arteryTAC(lp)<0
        arteryTAC(lp)=0;
    end
    addpoints(arteryProfile,PETframingMidPoint(lp),arteryTAC(lp));
    addpoints(cranialMixedZoneProfile,PETframingMidPoint(lp),cranialTAC(lp));
    addpoints(caudalMixedZoneProfile,PETframingMidPoint(lp),caudalTAC(lp));
    arteryProfile.DisplayName='Artery activity';
    brainProfile.DisplayName='Brain activity';
    cranialMixedZoneProfile.DisplayName='Cranial mixed zone activity';
    caudalMixedZoneProfile.DisplayName='Caudal mixed zone activity';
    legend(arteryProfile.DisplayName,brainProfile.DisplayName,cranialMixedZoneProfile.DisplayName,caudalMixedZoneProfile.DisplayName,'FontSize',20);
    drawnow;
end
arteryTAC(arteryTAC<0)=0;

interpArteryTAC=interp1(PETframingMidPoint,arteryTAC,1:3600,'pchip');
pvcOutputs.arteryTACarea=trapz(interpArteryTAC);
% Logic for performing partial volume correction based on activity distribution as a function of time.


% Only three activity curves should be considered and are crucial according
% to the current aspect of modelling 
framesContainingProminentSpillOutFromArtery=(round(arteryTAC,2)>round(brainTAC,2));% & (round(caudalTAC,2)>round(cranialTAC,2));
framesContainingProminentSpillOutFromBrain=(round(arteryTAC,2)<round(brainTAC,2));% & (round(caudalTAC,2)<round(cranialTAC,2)) & (round(arteryTAC,2)<round(cranialTAC,2));
% startFrame=find(framesContainingProminentSpillOutFromArtery);
% endFrame=find(framesContainingProminentSpillOutFromBrain);
% framesContainingEqualSpillOutFromBrainAndArtery=(startFrame(end)+1):(endFrame(end));
earlyFrames=find(framesContainingProminentSpillOutFromArtery);
lateFrames=find(framesContainingProminentSpillOutFromBrain);
% middleFrames=framesContainingEqualSpillOutFromBrainAndArtery;
% middleFrames=[earlyFrames(end)+1:lateFrames(1)-1];
plot(PETframingMidPoint(framesContainingProminentSpillOutFromArtery),arteryTAC(framesContainingProminentSpillOutFromArtery),'mo--','LineWidth',2)
hold on
plot(PETframingMidPoint(framesContainingProminentSpillOutFromBrain),arteryTAC(framesContainingProminentSpillOutFromBrain),'ko--','LineWidth',2)

%%%%%
% Perform correction for early frames.
% Logic behind this correction: Early frames have very high activity in the
% artery and very less activity in the brain. Therefore it is sufficient to
% perform spill-out correction only.
PVCProfile=animatedline('Color','cyan','Marker','o','LineWidth',3);
PVCProfile_1=animatedline('Color','blue','Marker','p','LineWidth',3);
PVCProfile_2=animatedline('Color','magenta','Marker','p','LineWidth',3);
frameSegmentOne=earlyFrames;
if isempty(frameSegmentOne)
    frameSegmentOne=1:26;
end
for lp=1:length(frameSegmentOne)
   indexOfInterest=frameSegmentOne(lp);
   spillOutCorrected=PET{indexOfInterest}./recoveryCoefficients{indexOfInterest};
    pvcTarget(indexOfInterest)=(mean(spillOutCorrected(targetMRmask{indexOfInterest}>0))./1000);
    [mask_1,mask_2]=VolumePartitioner(targetMRmask{indexOfInterest}>0);
    pvcTarget1(indexOfInterest)=(mean(spillOutCorrected(mask_1>0))./1000);
    pvcTarget2(indexOfInterest)=(mean(spillOutCorrected(mask_2>0))./1000);
    if pvcTarget(indexOfInterest)<0
        pvcTarget(indexOfInterest)=0;
        pvcTarget1(indexOfInterest)=0;
        pvcTarget2(indexOfInterest)=0;
    end
    addpoints(PVCProfile,PETframingMidPoint(indexOfInterest),pvcTarget(indexOfInterest));
    drawnow
    text(PETframingMidPoint(lp),pvcTarget(lp)+2,num2str(lp));

end

if PVC.maskSwitch=='c' 
    if isempty(frameSegmentOne)
        frameSegmentOne=1:26;
    else
    end
    frameSegmentOne(end)=26;
end

%% Partial volume correction for middle frames: Artery activity ~= Brain activity

for lp=(frameSegmentOne(end)):33
%for lp=(1):33
    fusedBrain=significantBrain{lp}.*PET{lp}; % sample the background brain activity from PET, which contributes to the spill-over in a significant manner.
    scaledBrain=performOtsu(fusedBrain,numberOfThresholds); % account for circumferential variability
    brainContribution=imgaussfilt3(scaledBrain,psfSigma); % calculate the spill-over factors.
    brainContribution=brainContribution.*~significantBrain{lp}; % isolate spill-over factors to mixed-zone, not to the brain.
    brainContribution=brainContribution.*~targetMRmask{lp}; % isolate spill-over factors to mixed-zone, not to the artery.
    modPET=PET{lp}; % Passing the original PET volume to a dummy PET variable.
    %relDiff=10;
    convergenceAchieved=false;
    [initialEstimate,newPET, newEstimate_1,newEstimate_2]=hardPVC(modPET,backgroundMantel{lp},numberOfThresholds,psfSigma,recoveryCoefficients{lp},double(targetMRmask{lp}));
    initialEstimate=initialEstimate*1000; % Initial crude estimate of the ICA, obtained from the previous step and conversion from kBq/mL to Bq/mL
    incrm=0;
    while convergenceAchieved==false
            incrm=incrm+1; 
        % clean mixed zone and recalculate
            modPET=PET{lp}-brainContribution;
            modPET(modPET<0)=0;
            scaledArtery=initialEstimate.*targetMRmask{lp};
            arteryContribution=imgaussfilt3(scaledArtery,psfSigma);
            arteryContribution=arteryContribution.*~targetMRmask{lp};
     %       arteryContribution=arteryContribution.*~significantBrain{lp};
            modPET=modPET-arteryContribution;
            modPET(modPET<0)=0;
            [newEstimate,newPET, newEstimate_1,newEstimate_2]=hardPVC(modPET,backgroundMantel{lp},numberOfThresholds,psfSigma,recoveryCoefficients{lp},double(targetMRmask{lp}));
            newEstimateVec(incrm)=newEstimate;
            relDiff=abs((((newEstimate*1000)-initialEstimate)./initialEstimate)*100);
            disp(['Processing Frame: ',num2str(lp),'! The relative difference is ',num2str(relDiff),'!'])
            initialEstimate=newEstimate*1000;
            relDiffVec(incrm)=relDiff;
            relDiffVec(relDiffVec==0)=[];
            if relDiff==0
               convergenceAchieved=true;
               disp(['Convergence achieved at iteration no:',num2str(incrm)])
               break
            end
            if length(relDiffVec)==3
               signChange=sign(diff(relDiffVec));
               minimaCheck=sum(signChange);
               if minimaCheck==0
                  convergenceAchieved = true;
                  disp(['Convergence achieved, minimum found at iteration no: ', num2str(incrm),' with a relative difference value of ',num2str(relDiffVec(2)),'!']);    
                  clear relDiffVec
               else
                  convergenceAchieved = false;
                  clear relDiffVec
                  relDiffVec(1)=relDiff; % move the last calculated relative difference to the first of
                  incrm=1;
               end
            end
    end   
        if newEstimate<0
            newEstimate=0;
        end
        if newEstimate_1<0
            newEstimate_1=0;
        end
        if newEstimate_2<0
            newEstimate_2=0;
        end
        pvcTarget(lp)=newEstimate;
        pvcTarget1(lp)=newEstimate_1;
        pvcTarget2(lp)=newEstimate_2;
        addpoints(PVCProfile,PETframingMidPoint(lp),pvcTarget(lp));
        addpoints(PVCProfile_1,PETframingMidPoint(lp),pvcTarget1(lp));
        addpoints(PVCProfile_2,PETframingMidPoint(lp),pvcTarget2(lp));
        drawnow
        text(PETframingMidPoint(lp),pvcTarget(lp)+2,num2str(lp));

end

%% Partial volume correction for middle frames: Brain activity >> Artery activity

for lp=34:length(PET)
    fusedBrain=significantBrain{lp}.*PET{lp};
    scaledBrain=performOtsu(fusedBrain,numberOfThresholds);
    brainContribution=imgaussfilt3(scaledBrain,psfSigma);
    brainContribution=brainContribution.*~significantBrain{lp};
    newBrainContribution=brainContribution;
    modPET=PET{lp};
    %relDiff=10;
    convergenceAchieved = false;
    [initialEstimate,newPET, newEstimate_1,newEstimate_2]=hardPVC(modPET,backgroundMantel{lp},numberOfThresholds,psfSigma,recoveryCoefficients{lp},double(targetMRmask{lp}));
    initialEstimate=initialEstimate*1000; % conversion from kBq/mL to Bq/mL
    incrm=0;
    while convergenceAchieved == false %relDiff>3 
        incrm=incrm+1; 
    % clean mixed zone and recalculate
        modPET=PET{lp}-newBrainContribution;
        modPET(modPET<0)=0;
        scaledArtery=initialEstimate.*targetMRmask{lp};
        arteryContribution=imgaussfilt3(scaledArtery,psfSigma);
        arteryContribution=arteryContribution.*~targetMRmask{lp};
        arteryContribution=arteryContribution.*~significantBrain{lp};
        modPET=modPET-arteryContribution;
        modPET(modPET<0)=0;
        [newEstimate,newPET, newEstimate_1,newEstimate_2]=hardPVC(modPET,backgroundMantel{lp},numberOfThresholds,psfSigma,recoveryCoefficients{lp},double(targetMRmask{lp}));
        newEstimateVec(incrm)=newEstimate;
        relDiff=abs((((newEstimate*1000)-initialEstimate)./initialEstimate)*100);
        disp(['Processing Frame: ',num2str(lp),'! The relative difference is ',num2str(relDiff),'!'])
        initialEstimate=newEstimate*1000;
        relDiffVec(incrm)=relDiff;
        relDiffVec(relDiffVec==0)=[]; 
        if relDiff==0
            convergenceAchieved=true;
            disp(['Convergence achieved at iteration no:',num2str(incrm)])
            break
        end
        if length(relDiffVec)==3
            signChange=sign(diff(relDiffVec));
            minimaCheck=sum(signChange);
            if minimaCheck==0
                convergenceAchieved = true;
                disp(['Convergence achieved, minimum found at iteration no: ', num2str(incrm),' with a relative difference value of ',num2str(relDiffVec(2)),'!']);
                clear relDiffVec
            else
                convergenceAchieved = false;
                clear relDiffVec
                relDiffVec(1)=relDiff; % move the last calculated relative difference to the first of  
                incrm=1;
            end
        end
    end   
    if newEstimate<0
        newEstimate=0;
    end
    if newEstimate_1<0
        newEstimate_1=0;
    end
    if newEstimate_2<0
        newEstimate_2=0;
    end
        %  close(gcf)
    pvcTarget(lp)=newEstimate;
    pvcTarget1(lp)=newEstimate_1;
    pvcTarget2(lp)=newEstimate_2;
    addpoints(PVCProfile,PETframingMidPoint(lp),pvcTarget(lp));
    addpoints(PVCProfile_1,PETframingMidPoint(lp),pvcTarget1(lp));
    addpoints(PVCProfile_2,PETframingMidPoint(lp),pvcTarget2(lp));
    drawnow
    text(PETframingMidPoint(lp),pvcTarget(lp)+2,num2str(lp));
    
end 
if isempty(pvcTarget1) 
   pvcTarget1=pvcTarget;
else
    if isnan(pvcTarget1)
        pvcTarget1=pvcTarget;
    else
    end
end
if isempty(pvcTarget2) 
    pvcTarget2=pvcTarget;
else
    if isnan(pvcTarget2)
        pvcTarget2=pvcTarget;
    else
    end
end
IDIFpostProcessingParameters.time=PETframingMidPoint./60;
IDIFpostProcessingParameters.P2Bratio=PVC.P2Bratio;
IDIFpostProcessingParameters.IDIF=pvcTarget;
[postProcessedIDIF]=IDIFpostProcessing(IDIFpostProcessingParameters);
IDIFpostProcessingParameters.IDIF=pvcTarget1;
[postProceesedIDIF_1]=IDIFpostProcessing(IDIFpostProcessingParameters);
IDIFpostProcessingParameters.IDIF=pvcTarget2;
[postProceesedIDIF_2]=IDIFpostProcessing(IDIFpostProcessingParameters);

time=1:3600;
plot(time,postProcessedIDIF,'b-','LineWidth',2)
%% // to uncomment
interpAIF=interp1((r.time_minutes_.*60),r.value_kBq_cc_,1:3600,'pchip');
interpAIF(isnan(interpAIF))=0;
[alignedIDIF,alignedAIF,delayPoints]=alignsignals((postProcessedIDIF),interpAIF);
[alignedIDIF1,alignedAIF,delayPoints]=alignsignals((postProceesedIDIF_1),interpAIF);
[alignedIDIF2,alignedAIF,delayPoints]=alignsignals((postProceesedIDIF_2),interpAIF);
figure,plot(time,alignedIDIF(time),'mo-');hold on
plot(time,alignedIDIF1(time),'ko-');
plot(time,alignedIDIF2(time),'yo-');
plot(time,alignedAIF(time),'bo-');
hLegend=legend('Image-derived input function','Mask-1: IDIF','Mask-2: IDIF','Arterial input function')
hLegend.FontSize=20;
xlabel('Time in seconds','FontSize',20);
ylabel('Activity in kBq/mL','FontSize',20);
aucAIF=trapz(time,alignedAIF(time));
aucIDIF=trapz(time,alignedIDIF(time));
rAUC=aucIDIF/aucAIF;
disp([' The similarity ratio is ',num2str(rAUC),'!'])
aucIDIF_1=trapz(time,alignedIDIF1(time));
rAUC_1=aucIDIF_1/aucAIF;
aucIDIF_2=trapz(time,alignedIDIF2(time));
rAUC_2=aucIDIF_2/aucAIF;
disp([' The similarity ratio is ',num2str(rAUC_1),'!'])
disp([' The similarity ratio is ',num2str(rAUC_2),'!'])

frameNumbers=1:length(PET); rotAlignedIDIF=alignedIDIF'; 
text(PETframingMidPoint,(rotAlignedIDIF(round(PETframingMidPoint))+2),num2str(frameNumbers'));

% // to uncomment 
%% Saving the analysis report. // to uncomment

title([PVC.subjectID,': The similarity ratio is ',num2str(rAUC), 'Mask_1: ',num2str(rAUC_1),' Mask_2: ',num2str(rAUC_2),'!'],'FontSize',20);
set(gcf,'position',get(0,'ScreenSize')) % expand figure to whole screen;
cd(PVC.whereToProcess);
saveas(gcf,[PVC.AC,PVC.subjectID,'IDIF_AIF','.png']);
IDIFtxtCreator(time,alignedIDIF(time),[PVC.AC,PVC.subjectID ,'_Aligned_IDIF']);
delayPoints=-delayPoints; % change the sign for matching the AIF to IDIF
modTime=(time+delayPoints);
IDIFtxtCreator(modTime,interpAIF,[PVC.AC,PVC.subjectID,'_DC_AIF']);
% // to uncomment
title([PVC.subjectID],'FontSize',20);
set(gcf,'position',get(0,'ScreenSize')) % expand figure to whole screen;
cd(PVC.whereToProcess);
IDIFtxtCreator(time,postProcessedIDIF,[PVC.AC,PVC.subjectID,'_native_IDIF']);
IDIFtxtCreator(time,postProceesedIDIF_1,[PVC.AC,PVC.subjectID,'_mask_1_IDIF']);
IDIFtxtCreator(time,postProceesedIDIF_2,[PVC.AC,PVC.subjectID,'_mask_2_IDIF']);
%IDIFtxtCreator(PETframingMidPoint,arteryTAC,[PVC.AC,'_No_MoCo_',PVC.subjectID,'_RAW_IDIF']);
%save('Artery_values',arteryTAC);
%% // to uncomment
pvcOutputs.IDIF=alignedIDIF(time);
pvcOutputs.Time=time;
pvcOutputs.rAUC1=rAUC_1;
pvcOutputs.rAUC2=rAUC_2;
pvcOutputs.rAUC=rAUC;
pvcOutputs.RMSE = sqrt((mean((alignedAIF(time)-alignedIDIF(time)).^2)));
pvcOutputs.APD=abs((aucIDIF-aucAIF)/aucAIF)*100;
disp(['Absolute percentage difference: ',num2str(pvcOutputs.APD)]);
pvcOutputs.rawAUC=(pvcOutputs.arteryTACarea./aucAIF)
disp(['RMSE: ', num2str(pvcOutputs.RMSE)]); 
% to uncomment
end

% function which performs the second pvc after the background is modelled accurately.
function [lateFramesPVC,newPET, newEstimate_1,newEstimate_2]=hardPVC(PETimg,backgroundRegion,numberOfThresholds,psfSigma,recoveryCoefficients,targetMRmask)
        fusedBG=PETimg.*backgroundRegion;
        scaledBackground=performOtsu(fusedBG,numberOfThresholds);
        spillInToArtery=imgaussfilt3(scaledBackground,psfSigma);
        spillInCorrectedPET=PETimg-spillInToArtery;
        spillInCorrectedPET(spillInCorrectedPET<0)=0;
        spillOutCorrectedPET=spillInCorrectedPET./recoveryCoefficients;
        spillOutCorrectedPET(isinf(spillOutCorrectedPET))=0;
        spillOutCorrectedPET(isnan(spillOutCorrectedPET))=0;
        newPET=spillOutCorrectedPET;
        pvcTarget=mean(spillOutCorrectedPET(targetMRmask>0))./1000;
        lateFramesPVC=pvcTarget;
        [mask1,mask2]=VolumePartitioner(targetMRmask>0);
        newEstimate_1=mean(spillOutCorrectedPET(mask1>0))./1000;
        newEstimate_2=mean(spillOutCorrectedPET(mask2>0))./1000;
end
