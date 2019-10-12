%-------------------------------------------------------------------------%
% This is a function which allocates the dynamic PET frames  to their
% respective MR navigators based on least temporal difference between their
% time stamps.
% 
% Author: Lalith Kumar Shiyam Sundar, M.Sc.
% Date: January 26, 2018
% Quantitative Imaging and Medical Physics
% 
% Inputs: AMP.pathOfDicomPET= path of the DICOM PET data.
%         AMP.pathOfDicomMRnavigators = path of the DICOM MR navigators.
%           
% 
% Usage:
%       MRnavigatorCorrespondenceToPET=allocateMRnavigatorsToPET(AMP)
%
%-------------------------------------------------------------------------%
%                           Program Start 
%-------------------------------------------------------------------------%

function MRnavigatorCorrespondenceToPET=allocateMRnavigatorsToPET(AMP)


[PETframingMidPoint]=getPETmidTime;


%% Obtain MR time stamps.

cd(AMP.pathOfDicomMRnavigators); % MR navigators will have sub-folders;
foldersInside=dir;
foldersInside=foldersInside(arrayfun(@(x) x.name(1),foldersInside) ~= '.');
for lp=1:length(foldersInside)
    prepFoldersForSort{lp,:}=foldersInside(lp).name;
end
sortedFolders=natsort(prepFoldersForSort);

for lp=1:length(sortedFolders)
    temp=[cd,filesep,sortedFolders{lp,:}];
    cd(temp);
    filesInside=dir; % read the folders
    filesInside=filesInside(arrayfun(@(x) x.name(1), filesInside) ~= '.'); % read the files in the path;
    MRdcmInfo=dicominfo(filesInside(1).name);
    MRframes(lp)=lp;
    MRdateStamps(lp)=datetime([MRdcmInfo.AcquisitionDate,MRdcmInfo.AcquisitionTime(1:6)],'InputFormat','yyyyMMddHHmmss');
    MRacquisitionTime(lp)=abs(seconds(MRdateStamps(lp)-MRdateStamps(1)));
    cd .. 
end
MRacquisitionTime=MRacquisitionTime';
MRframes=MRframes';
MRdateStamps=MRdateStamps';
MRtimingInfo=table(MRframes,MRdateStamps);


%% Calculate elapsed time difference between the PET frames and each MR navigator
figure('units','normalized','outerposition',[0 0 1 1]);

for lp=1:length(MRframes)
    timeDifference(lp,:)=abs((PETframingMidPoint-MRacquisitionTime(lp)));
    plot(timeDifference(lp,:),'o-','LineWidth',2);
    hold on
    pause(0.1)
end
xlabel('PET frames','FontSize',20)
ylabel('Time difference between the MR navigators and PET frames','FontSize',20)
for lp=1:length(PETframingMidPoint)
    temporalMinimum=find(timeDifference(:,lp)==min(timeDifference(:,lp)));
    PETlinkedNavigator(lp)=temporalMinimum(1);
    selectedMinima=min(timeDifference(:,lp));
    plot(lp,selectedMinima,'ro','LineWidth',7,'MarkerFaceColor',[0.5 0.5 0.5]);
    text(lp+0.25,selectedMinima+50,num2str(PETlinkedNavigator(lp)),'FontSize',14);
    pause(0.1)
end
    
MRnavigatorCorrespondenceToPET(1,:)=1:length(PETframingMidPoint);
MRnavigatorCorrespondenceToPET(2,:)=PETlinkedNavigator;
MRnavigatorCorrespondenceToPET=MRnavigatorCorrespondenceToPET';
h=msgbox('MR navigators linked with PET frames based on least temporal difference!','Success');

% Write a text file which contains the MR navigator link to the PET

fid=fopen('PET_correspondence.txt','w');
fprintf(fid,'%d \r\n',[MRnavigatorCorrespondenceToPET(:,1)]);
fclose(fid);
fid=fopen('MR_correspondence.txt','w');
fprintf(fid,'%d \r\n',[MRnavigatorCorrespondenceToPET(:,2)]);
fclose(fid);

end


