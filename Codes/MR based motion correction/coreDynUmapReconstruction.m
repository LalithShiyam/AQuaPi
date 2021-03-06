%% This function is specifically designed for the IDIF protocol in AKH/QIM. This wrapper will rotate the attenuation maps based on the transformation matrices obtained from the MR Navigators. 
%% and perform the dynamic reconstructions using the JSrecon e7 toolkit.
%   
%  Author: Lalith Kumar Shiyam Sundar, M.Sc.
%  Date: Mar 23rd, 2018
%  Quantitative Imaging and Medical Physics, Medical University of Vienna.
%
%  Inputs: 
%       DURinputs.subjectID= ID of the subject.
%       DURinputs.pathOfStaticACmaps= path of the static u-maps (can be
%       more than 1)
%       DURinputs.pathOfReconParamFile= path of the reconstruction
%       parameter file
%       DURinputs.pathOfDicomMRnavigators= path where the DICOM MR
%       navigators are stored
%       DURinputs.pathToListModePET= path pointing towards the PET
%       list-mode data.
%       DURinputs.pathToCopyReconPETdata = path to copy the reconstructed
%       pet data.
%       
%  Outputs: 
%       Dynamic reconstructions of the PET data with motion induced AC
%       maps.
%%------------------------------------------------------------------------%
%                               Program start
%%------------------------------------------------------------------------%

function []=coreDynUmapReconstruction(DURinputs)

%% Deriving motion imposed attenuation maps.

% Obtaining the motion imposed attenuation maps using the MR navigators,
% the output u-maps are in DICOM and are JSRecon 12 compatible (version
% - xx)

% passing the global inputs as local function inputs - purely for the sake
% of readability (*for new programmers)


cd(DURinputs.pathOfStaticACmaps);
genMIUmaps.pathOfStaticAC=DURinputs.pathOfStaticACmaps;
genMIUmaps.pathOfMRnavigators=DURinputs.pathOfDicomMRnavigators;
cd(genMIUmaps.pathOfMRnavigators); % go to the MR navigator folder 
fNames=dir; % read the folders inside the mr navigator folder
fNames=fNames(arrayfun(@(x) x.name(1), fNames) ~= '.'); % read volunteer folders
if isempty(fNames) % check if the folder has mr navigators % if there are none, simply proceed.
   disp('No MR navigators found for this subject!'); 
else 
   disp('MR navigators found, rotating attenuation maps...'); 
   generateMotionImposedUmaps(genMIUmaps); % if there are mr navigators, use the idifRotateACmaps program to obtain motion induced ac maps
end

%% reading the paths of dynamic u-maps

cd(genMIUmaps.pathOfStaticAC); 
cd .. 
dynamicFolder=dir('Dynamic*');
pathOfDynamicFolder=[pwd,filesep,dynamicFolder.name];
cd(pathOfDynamicFolder)
dynamicFolders=dir; dynamicFolders=dynamicFolders(arrayfun(@(x) x.name(1), dynamicFolders)  ~= '.');
for lp=1:length(dynamicFolders)
    pathToDynUmapList{lp}=[pwd,filesep,dynamicFolders(lp).name]; % storing the path of the dynamic folders containing the rotated attenutaion maps.
end


%% Linking the prospective PET frames mid-point with the MR navigators 

AMP.pathOfReconParam=DURinputs.pathOfReconParamFile; % refer to the allocateMRnavigatorsToPET program documentation 
AMP.pathOfDicomMRnavigators=DURinputs.pathOfDicomMRnavigators; % refer to the allocateMRnavigatorsToPET program documentation 
AMP.pathToStoreTxt=AMP.pathOfReconParam; % where do you want to store the linking text file.
MRnavigatorCorrespondenceToPET=allocateMRnavigatorsToPET(AMP); % program which links prospective PET frame mid-point with the MR navigators.

%% Reconstructing the PET list-mode based on dynamic u-maps.

subjectID = DURinputs.subjectID; % getting the subject ID for the powershell scripts.
subjectID = regexprep(subjectID, ' ', '` ');
subjectID=['"',subjectID,'"'];
processingFolder='Desktop'; % this should always be desktop (purely for sorting purposes in windows - don't ask windows is stupid).
processingFolder=['"',processingFolder,'"'];
pathOfListModePET=DURinputs.pathToListModePET; % points to the path of the list mode data.
pathOfListModePET= regexprep(pathOfListModePET, ' ', '` ');
pathOfListModePET=['"',pathOfListModePET,'"'];
cd(DURinputs.pathOfReconParamFile); % goes to the path of recon param file
paramFile=dir('*s.txt'); %looks for the reconstruction parameter file 
pathOfParamFile=[DURinputs.pathOfReconParamFile,filesep,paramFile.name];
pathOfParamFile=regexprep (pathOfParamFile,' ','` ');
pathOfParamFile=['"',pathOfParamFile,'"'];
mrCorrespondenceLink=dir('*pondence.txt'); % looks for the pet-mr correspondence link 
pathOfPETMRlink=[DURinputs.pathOfReconParamFile,filesep,mrCorrespondenceLink.name];
pathOfPETMRlink=regexprep(pathOfPETMRlink,' ','` ');
pathOfPETMRlink=['"',pathOfPETMRlink,'"'];
pathToCopyReconPETdata=DURinputs.pathToCopyReconPETdata; % where to copy the recon pet data.
pathToCopyReconPETdata=regexprep(pathToCopyReconPETdata,' ','` ');
pathToCopyReconPETdata=['"',pathToCopyReconPETdata,'"'];
for lp=1:length(dynamicFolders)
    attenuationMapToCopy=dynamicFolders(lp).name; % variable which contains the type of attenuation map.
    attenuationMapToCopy=['"',attenuationMapToCopy,'"'];
    pathToDynUmap=pathToDynUmapList{lp};
    pathToDynUmap=regexprep(pathToDynUmap,' ','` ');
    pathToDynUmap=['"',pathToDynUmap,'"'];
    % This "if" construct checks for indication, if its the same patient  with different ac maps
    if lp==1
       samePatient=0; % indicates its a new patient
    else
       samePatient=1; % indicates its the same patient
    end
    samePatient=num2str(samePatient);
    samePatient=['"',samePatient,'"'];
   % Calling the reconstruction powershell script for volunteers
    % without MR navigators.
    cd(genMIUmaps.pathOfMRnavigators)
    fNames=dir;
    fNames=fNames(arrayfun(@(x) x.name(1), fNames) ~= '.'); % read volunteer folders
    prefixPowShellScript=['!powershell -executionPolicy bypass ". ']; % with dot sourcing.
    if isempty(fNames)
       staticPS=which('StaticRecon.ps1');
       pathOfStaticPS=[staticPS,'; '];
       pathOfStaticPS=regexprep(pathOfStaticPS,' ','` ');
       inputArgString=['StaticRecon ',subjectID,' ',processingFolder,' ',attenuationMapToCopy,' ',pathToDynUmap,' ',pathOfListModePET,' ',pathOfParamFile,' ',pathOfPETMRlink,' ',pathToCopyReconPETdata,' ',samePatient]; 
       stringToRun=[prefixPowShellScript,pathOfStaticPS,inputArgString];
       eval(stringToRun);
    else
       dynamicPS=which('MoCoRecon.ps1');
       pathOfDynamicPS=[dynamicPS,'; '];
       pathOfDynamicPS=regexprep(pathOfDynamicPS,' ','` ');
       inputArgString=['MoCoRecon ',subjectID,' ',processingFolder,' ',attenuationMapToCopy,' ',pathToDynUmap,' ',pathOfListModePET,' ',pathOfParamFile,' ',pathOfPETMRlink,' ',pathToCopyReconPETdata,' ',samePatient]; 
       stringToRun=[prefixPowShellScript,pathOfDynamicPS,inputArgString];
       eval(stringToRun); 
    end
end
cd(pathToCopyReconPETdata);
cd ..
pathOfNifti=pwd;
constructDynPETseries(pathToCopyReconPETdata,pathOfNifti);
end



    
        
    