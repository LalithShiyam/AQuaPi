
%% Input parameters:
    PVC.pathOfPET='';
    PVC.pathOfhighResMask='';
    PVC.pathOfBrainMasks= '';
    PVC.pathOfpMoCohighResMask='';
    PVC.pathOfpMoCoBrainMasks= '';
    PVC.pathOfReconParam= '';
    PVC.pathToMRcorrespondenceFile= '';
    
%% //to uncomment
    [dateVector]=DCMinfoExtractor(subjectID);

    tLower1=datetime(2016,04,11); tUpper1=datetime(2016,06,17)
    tLower2=datetime(2016,06,25); tUpper2=datetime(2016,10,13)
    tLower3=datetime(2016,10,24); tUpper3=datetime(2017,1,12)
    tLower4=datetime(2017,01,23); tUpper4=datetime(2017,04,25)

    if isbetween(dateVector(lp,:),tLower1,tUpper1)
       PVC.crossCalibrationFactor=1.141
    end
    if isbetween(dateVector(lp,:),tLower2,tUpper2)
       PVC.crossCalibrationFactor=1.126;
    end
    if isbetween(dateVector(lp,:),tLower3,tUpper3)
       PVC.crossCalibrationFactor=1.111;
    end
    if isbetween(dateVector(lp,:),tLower4,tUpper4)
       PVC.crossCalibrationFactor= 1.0880
    end
%% // to uncomment

   %PVC.crossCalibrationFactor=1.115;
    % calculate plasma to blood ratio for individual subject
    volunteerID=shortSubjectID{lp}
%     %%// to uncomment
    pathOfAIF='/Users/lalithsimac/Documents/Input functions/Arterial input functions';% arterial input function
    cd(pathOfAIF);
    AIFtable=dir([volunteerID,'*.txt']);
    r=readtable(AIFtable.name);
%// to uncomment

    %% // to uncomment

    pathOfBIF='/Users/lalithsimac/Documents/Input functions/Blood input functions'; % whole blood input function
    cd(pathOfBIF);
    BIFtable=dir([volunteerID,'*.txt'])
    q=readtable(BIFtable.name);
    P2Bratio=r.value_kBq_cc_./q.whole_blood_kBq_cc_; %plasma/wholeblood ratio
    P2Bratio(isnan(P2Bratio))=[];
    P2Bratio=mean(P2Bratio);
    PVC.P2Bratio=P2Bratio;
% % // to  uncomment
 PVC.P2Bratio=1.06

    PVC.psfFWHM= 5

    PVC.petMoCo=false;
    PVC.AC=[acMaps{ilp}]%,'_inferior'];
    PVC.whereToProcess= [pathSuffix,subjectID{lp},'/Processed data'];
    % to comment out
   % r=[]; 
    % to comment out
    [pvcOutputs] = tempPVC2(PVC,r);
%     if strcmp(acMaps,'RESOLUTE')
%         resoluteAUC1(lp)=pvcOutputs.rAUC1;
%         resoluteAUC2(lp)=pvcOutputs.rAUC2;
%         resoluteAUC(lp)=pvcOutputs.rAUC;
%     else
%         if strcmp(acMaps,'PCT')
%             PCTauc1(lp)=pvcOutputs.rAUC1;
%             PCTauc2(lp)=pvcOutputs.rAUC2;
%             PCTauc(lp)=pvcOutputs.rAUC;
%         end
%     end
    %rawArtery{lp}(ilp)=pvcOutputs.arteryTACarea;
%     rawAUC{lp}(ilp)= pvcOutputs.rawAUC;
end
end
% resoluteAUC1=resoluteAUC1';
% resoluteAUC2=resoluteAUC2';
% resoluteAUC=resoluteAUC';
% PCTauc=PCTauc';
% PCTauc1=PCTauc1';
% PCTauc2=PCTauc2';
% T=table(subjectID,resoluteAUC,resoluteAUC1,resoluteAUC1,resoluteAUC2,PCTauc,PCTauc1,PCTauc2);
% cd(whereToStoreResults)
% writetable(T);
toc
