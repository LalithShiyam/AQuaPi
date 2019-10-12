%% This function creates motion imposed u-maps obtained from the motion 
% vectors derived from the MR navigators.

% Input: 
%           genMIUmaps.pathOfStaticAC= path of the static u-maps.
%           genMIUmaps.pathOfMRnavigators = path of the mr navigators. 
%           

%% Program start.

% create the folder structure 
function[]=generateMotionImposedUmaps(genMIUmaps)
cd(genMIUmaps.pathOfStaticAC)
cd .. 
whereAmI=cd
cd(whereAmI);
mkdir('Dynamic AC')
mkdir('Working AC')
pathOfDynAC=[whereAmI,filesep,'Dynamic AC'];
pathOfWorkingArea=[whereAmI,filesep,'Working AC'];
%% Read the different type of attenuation maps
cd(genMIUmaps.pathOfStaticAC);
ACnames=dir;
ACnames=ACnames(arrayfun(@(x) x.name(1), ACnames) ~= '.'); % read volunteer folders
parfor lp=1:length(ACnames)
    cd(pathOfDynAC)
    mkdir(ACnames(lp).name);
    cd(pathOfWorkingArea)
    mkdir(ACnames(lp).name);
    disp([ACnames(lp).name,' Created in',pathOfDynAC,' and ',pathOfWorkingArea]);
end

%% calculate the number of mr navigators and create the same number of dynamic u-maps.
cd(genMIUmaps.pathOfMRnavigators);
mrNavFiles=dir;
mrNavFiles=mrNavFiles(arrayfun(@(x) x.name(1), mrNavFiles) ~= '.'); % read volunteer folders
numberOfMRnav=length(mrNavFiles);
% Copy the static folders to dynamic folders based on the number of mr
% navigators present
cd(pathOfDynAC);
wb=timebar('Preparing the buffer u-maps...');
for olp=1:length(ACnames)
    parfor ilp=1:(numberOfMRnav)
        staticPath=[genMIUmaps.pathOfStaticAC,filesep,ACnames(olp).name];
        if strcmp(ACnames(olp).name,'RESOLUTE')
            staticPath=[genMIUmaps.pathOfStaticAC,filesep,'UTE'];
        else
        end
        ACseriesName=[ACnames(olp).name,'_uMap_',num2str(ilp)];
        folderOfInterest= dir([ACnames(olp).name,'*']);
        destinationPath=[pathOfDynAC,filesep,folderOfInterest.name,filesep,ACseriesName];
        copyfile (staticPath,destinationPath)
    end
    timebar(wb,olp/length(ACnames));
end
close(wb)

%% convert the static dicom AC files to nifti 

parfor lp=1:length(ACnames)
    cd([genMIUmaps.pathOfStaticAC,filesep,ACnames(lp).name])
    convertDicomtoNii(cd,cd)
    oldFile=dir('*.nii');
    newFile=[ACnames(lp).name,'.nii'];
    movefile(oldFile.name,newFile);
    movefile(newFile,[pathOfWorkingArea,filesep,ACnames(lp).name]);
    disp(['Converted DICOM files in ',cd,' to Nifti']);
end

%% Apriori knowledge: BD, CT, DIXON, PCT are in the same space. 

% Store all teh static nifty ac masks in the mask cell fro coregistration
incr=0;
for lp=1:length(ACnames)
    if strcmp(ACnames(lp).name,'BD') || strcmp(ACnames(lp).name,'CT') || strcmp(ACnames(lp).name,'DIXON') || strcmp(ACnames(lp).name,'PCT')
        incr=incr+1;
        CoregInputs.MaskImgPath{incr,:}=[pathOfWorkingArea,filesep,ACnames(lp).name,filesep,[ACnames(lp).name,'.nii']];
    else
    end
end

%% convert the mr navigators from their space to the dixon space. 

for lp=1:length(ACnames)
    if strcmp(ACnames(lp).name,'CT') || strcmp(ACnames(lp).name,'DIXON') || strcmp(ACnames(lp).name,'PCT')
        disp(['Reading files in ',genMIUmaps.pathOfStaticAC,filesep,ACnames(lp).name])
        cd([genMIUmaps.pathOfStaticAC,filesep,ACnames(lp).name])
        dixonSpaceUmap=readImages(cd);
        break
    else
    end
end
cd(pathOfWorkingArea)
mkdir('Dixon-space-navigators');
pathOfDixNav=[pathOfWorkingArea,filesep,'Dixon-space-navigators'];
cd(genMIUmaps.pathOfMRnavigators)
for lp=1:length(mrNavFiles)
    mrNavFilesToSort{lp,:}=mrNavFiles(lp).name;
end
sortedmrNav=natsort(mrNavFilesToSort);
parfor lp=1:length(sortedmrNav)
    cd([genMIUmaps.pathOfMRnavigators,filesep,sortedmrNav{lp}]);
    mrNav=readImages(cd)
    outMRnav=transformVoxelSpace(mrNav,dixonSpaceUmap);
    cd(pathOfDixNav)
    writeDicomImage(outMRnav,['Transform_MRnav_Dix_',num2str(lp)]);
    cd([pathOfDixNav,filesep,['Transform_MRnav_Dix_',num2str(lp)]]);
    convertDicomtoNii(cd,cd);
    oldFile=dir('*.nii');
    movefile(oldFile.name,['Transform_MRnav_Dix_',num2str(lp),'.nii'])
    movefile(['Transform_MRnav_Dix_',num2str(lp),'.nii'],pathOfDixNav);
    disp(['Transform_MRnav_Dix_',num2str(lp),'.nii',' created!']);    
end


%% perform the coregistration between mr navigators and apply them to the mask;
cd([pathOfDixNav])
tranFiles=dir('Transform_MRnav_Dix_*nii');
for lp=1:length(tranFiles)
    tranFilesToSort{lp,:}=tranFiles(lp).name;
end
sortedTF=natsort(tranFilesToSort)
CoregInputs.SourceImgPath=sortedTF{1};
CoregInputs.Interp=1; % trilinear interpolation 
for lp=2:numberOfMRnav
    disp(['Registering ',sortedTF{1}, ' to', sortedTF{lp}]);
    CoregInputs.Prefix=['First_To_',num2str(lp),'_'];
    CoregInputs.RefImgPath=[pathOfDixNav,filesep,sortedTF{lp}];
    Coregistration_job(CoregInputs)
end


%%
%% convert the mr navigators from their space to the dixon space. 

% for lp=1:length(ACnames)
%     if strcmp(ACnames(lp).name,'UTE') 
%         disp(['Reading files in ',genMIUmaps.pathOfStaticAC,filesep,ACnames(lp).name])
%         cd([genMIUmaps.pathOfStaticAC,filesep,ACnames(lp).name])
%         uteSpaceUmap=readImages(cd);
%         break
%     else
%     end
% end
% cd(pathOfWorkingArea)
% mkdir('ute-space-navigators');
% pathOfUteNav=[pathOfWorkingArea,filesep,'ute-space-navigators'];
% cd(genMIUmaps.pathOfMRnavigators)
% for lp=1:length(mrNavFiles)
%     mrNavFilesToSort{lp,:}=mrNavFiles(lp).name;
% end
% sortedmrNav=natsort(mrNavFilesToSort);
% parfor lp=1:length(sortedmrNav)
%     cd([genMIUmaps.pathOfMRnavigators,filesep,sortedmrNav{lp}]);
%     mrNav=readImages(cd)
%     outMRnav=transformVoxelSpace(mrNav,uteSpaceUmap);
%     cd(pathOfUteNav)
%     writeDicomImage(outMRnav,['Transform_MRnav_UTE_',num2str(lp)]);
%     cd([pathOfUteNav,filesep,['Transform_MRnav_UTE_',num2str(lp)]]);
%     convertDicomtoNii(cd,cd);
%     oldFile=dir('*.nii');
%     movefile(oldFile.name,['Transform_MRnav_UTE_',num2str(lp),'.nii'])
%     movefile(['Transform_MRnav_UTE_',num2str(lp),'.nii'],pathOfUteNav);
%     disp(['Transform_MRnav_UTE_',num2str(lp),'.nii',' created!']);    
% end
% clear CoregInputs % processing for ute
% incr=0;
% for lp=1:length(ACnames)
%     if strcmp(ACnames(lp).name,'UTE') || strcmp(ACnames(lp).name,'RESOLUTE') 
%         incr=incr+1;
%         CoregInputs.MaskImgPath{incr,:}=[pathOfWorkingArea,filesep,ACnames(lp).name,filesep,[ACnames(lp).name,'.nii']];
%     else
%     end
% end
% %% perform the coregistration between mr navigators and apply them to the mask;
% cd([pathOfUteNav])
% tranFiles=dir('Transform_MRnav_UTE_*nii');
% for lp=1:length(tranFiles)
%     tranFilesToSort{lp,:}=tranFiles(lp).name;
% end
% sortedTF=natsort(tranFilesToSort)
% CoregInputs.SourceImgPath=sortedTF{1};
% CoregInputs.Interp=1; % trilinear interpolation 
% for lp=2:numberOfMRnav
%     disp(['Registering ',sortedTF{1}, ' to', sortedTF{lp}]);
%     CoregInputs.Prefix=['First_To_',num2str(lp),'_'];
%     CoregInputs.RefImgPath=[pathOfUteNav,filesep,sortedTF{lp}];
%     Coregistration_job(CoregInputs)
% end


%% Push binary data into their respective dicom 
wb=timebar('Pushing rotated binary data to their respective DICOM AC maps...');
for olp=1:length(ACnames)
    if strcmp(ACnames(olp).name,'RESOLUTE')
        cd([pathOfWorkingArea,filesep,ACnames(olp).name])
        niftyVolume=spm_read_vols(spm_vol([ACnames(olp).name,'.nii']));
        searchString=(['*','_','1']);
        cd([pathOfDynAC,filesep,ACnames(olp).name]);
        dicomFolder=dir(searchString);
        pathOfDicom=[pathOfDynAC,filesep,ACnames(olp).name,filesep,dicomFolder.name];
        pushDataToDicom(pathOfDicom,niftyVolume);
        disp(['Pushing ',[ACnames(olp).name,'.nii'],' ',pathOfDicom,'!']);
    end
    for ilp=2:numberOfMRnav
        cd([pathOfDynAC,filesep,ACnames(olp).name])
        searchString=([ACnames(olp).name,'*','_',num2str(ilp)]);
        dicomFolder=dir(searchString);
        pathOfDicom=[pathOfDynAC,filesep,ACnames(olp).name,filesep,dicomFolder.name];
        cd([pathOfWorkingArea,filesep,ACnames(olp).name]);
        searchString=(['*','_',num2str(ilp),'_*nii']);
        niftyFile=dir(searchString);
        niftyVolume=spm_read_vols(spm_vol(niftyFile.name));
        pushDataToDicom(pathOfDicom,niftyVolume);
        disp(['Pushing ',niftyFile.name,' ',pathOfDicom,'!']);
    end
    timebar(wb,olp/length(ACnames));
end
close(wb)
end
%% Copy first resolute nifty volume to resolute-ute folder 

        
    
