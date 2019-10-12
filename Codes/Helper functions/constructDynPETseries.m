% This function constructs a nifti dynamic pet series from the split dicom
% pet images. The MoCoRecon.ps1 spits out chunked versions of the dynamic
% pet series. Since i am lazy to construct a dynamic pet series in dicom by
% changing the dicom tags or doing something fancy. I am going to convert
% everything to nifty and just rename their names. I know, i am lazy.
%
% Author: Lalith kumar shiyam sundar, M.Sc.
% Date: 27.04.2018
% 
% Inputs: pathOfSplitDicomPETSeries : Point the path where the split dicom
% series is present
% Output: Folder "Nifti Dynamic PET series"  
%
%------------------------------------------------------------------------%
%                           PROGRAM START
%------------------------------------------------------------------------%

% go the folder which contains the split dicom pet series 

function []=constructDynPETseries(pathOfSplitDicomPETSeries,whereToStore)
cd(pathOfSplitDicomPETSeries);
[name,~]=strsplit(pathOfSplitDicomPETSeries,filesep);
name=char(name(end));
niftiFolder=['Nifti_',name];
% read in the sub folders inside the split dicom series 

splitFolders=dir;
splitFolders=splitFolders(arrayfun(@(x) x.name(1), splitFolders) ~= '.'); 

% loop into a string container to perform natural sorting because matlab
% sorting sucks sometime

parfor lp=1:length(splitFolders)
    folderNames{lp,:}=splitFolders(lp).name;
end

sortedFolderName=natsort(folderNames); % sorting the folder names

cd(whereToStore)
mkdir(niftiFolder);
whereToStoreNiftiSeries=[pwd,filesep,niftiFolder];
%%
% conversion of the dicom pet images to nifty and renaming them based on a
% series convention.
cd(pathOfSplitDicomPETSeries)
frameCounter=0; 
incr=0;
for olp=1:length(sortedFolderName)
    tempPath=[pwd,filesep,sortedFolderName{olp,:}];
    [petFrames,~]=DynPETload(tempPath);
    for ilp=1:length(petFrames)
        incr=incr+1;
        dicomPETMeanValue(incr)=mean(petFrames{ilp}(:));
    end
    cd(tempPath);
    convertDicomtoNii(cd,cd);
    niftiFiles=dir('*nii');
    niftiFiles=niftiFiles(arrayfun(@(x) x.name(1), niftiFiles) ~= '.');
    for ilp=1:length(niftiFiles)
        niftiFilesStrings{ilp,:}=niftiFiles(ilp).name;
    end
    sortedNiftiFiles=natsort(niftiFilesStrings); % sorting the nifti files
    for ilp=1:length(sortedNiftiFiles)
        frameCounter=frameCounter+1;
        movefile(sortedNiftiFiles{ilp,:},['PET_Frame_',num2str(frameCounter),'.nii']);
    end
    clear niftiFilesStrings sortedNiftiFiles niftiFiles
    movefile('*.nii',whereToStoreNiftiSeries)

    cd .. 
end


%% checking if the converted nifti images are representative of the dynamic series
% check if the nifti images reflects the tracer dynamics

cd(whereToStoreNiftiSeries)
niftiFiles=dir('*nii');
niftiFiles=niftiFiles(arrayfun(@(x) x.name(1), niftiFiles) ~= '.');
parfor lp=1:length(niftiFiles)
    niftiFilesStrings{lp,:}=niftiFiles(lp).name;
end
sortedNiftiFiles=natsort(niftiFilesStrings);
figure,
sub1=subplot(2,1,1);
%axis([0 21 -1.1 1.1]);
h1=animatedline('Color','blue','Marker','o','LineWidth',3);
sub2=subplot(2,1,2)
%axis([0 21 -1.1 1.1]); 
h2=animatedline('Color','black','Marker','o','LineWidth',3);
for lp=1:length(sortedNiftiFiles)
    imgVol=spm_read_vols(spm_vol(sortedNiftiFiles{lp}));
    frameMeanValue(lp)=mean(imgVol(:));
    addpoints(h1,lp,frameMeanValue(lp));
    set(get(sub1,'XLabel'),'String','Nifti frames')
    set(get(sub1,'YLabel'),'String','Intensity');
    addpoints(h2,lp,dicomPETMeanValue(lp));
    set(get(sub2,'XLabel'),'String','DICOM frames')
    set(get(sub2,'YLabel'),'String','Intensity');
    drawnow;
end
%%
a=axes;
testOfEquality=mean(round(frameMeanValue-dicomPETMeanValue));
if testOfEquality==0
    t1=title('Perfect match between DICOM and NIFTI series');
    a.Visible='off';
    t1.Visible='on';
else
    t1=title('No match between DICOM and NIFTI series');
    a.Visible='off';
    t1.Visible='on';
end
end
