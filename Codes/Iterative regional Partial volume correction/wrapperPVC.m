%% Purpose
% This is a wrapper script for running the 'coreIterativePVC.m' program, 
% which is responsible for recovering the apparent activity in the carotid
% in an iterative manner.
%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna
%
%% Input parameters:
%     PVC.pathOfPET=path of reconstructed PET;
%     PVC.pathOfhighResMask=path of high resolution TOF-MRA segmented mask;
%     PVC.pathOfBrainMasks= path of the T1-MR brain mask;
%     PVC.pathOfpMoCohighResMask=path of the PET moco mr masks ;
%     PVC.pathOfpMoCoBrainMasks= path of the pet moco brain masks;
%     PVC.pathOfReconParam= path of the reconstruction parameters;
%     PVC.pathToMRcorrespondenceFile= path of the mr correspondence file generated from MR MoCO;
%     PVC.whereToProcess= path where to process the PVC;

%% Input parameters:
    PVC.pathOfPET='';
    PVC.pathOfhighResMask='';
    PVC.pathOfBrainMasks= '';
    PVC.pathOfpMoCohighResMask='';
    PVC.pathOfpMoCoBrainMasks= '';
    PVC.pathOfReconParam= '';
    PVC.pathToMRcorrespondenceFile= '';
    PVC.whereToProcess='';
%% Post-processing parameters, which needed to be updated based on experiments.

    [dateVector]=getDatesFromDcmInfo(PVC.pathOfPET)

    % Automatically calibrate the crosscalibration factor based on date.
    
    tLower1=datetime(2016,04,11); tUpper1=datetime(2016,06,17)
    tLower2=datetime(2016,06,25); tUpper2=datetime(2016,10,13)
    tLower3=datetime(2016,10,24); tUpper3=datetime(2017,1,12)
    tLower4=datetime(2017,01,23); tUpper4=datetime(2017,04,25)
    
    %---------------------------------------------------------------------%
    % For every new block of cross-calibration interval, the bottom code
    % needs to be appended in the following syntax:
    
    % if isbetween(datevector,tLowerX,tUpperX)
    %  PVC.crossCalibrationFactor=X.XXX
    % end
    %---------------------------------------------------------------------%
    
    if isbetween(dateVector,tLower1,tUpper1)
       PVC.crossCalibrationFactor=1.141
    end
    if isbetween(dateVector,tLower2,tUpper2)
       PVC.crossCalibrationFactor=1.126;
    end
    if isbetween(dateVector,tLower3,tUpper3)
       PVC.crossCalibrationFactor=1.111;
    end
    if isbetween(dateVector,tLower4,tUpper4)
       PVC.crossCalibrationFactor= 1.0880
    end

%% Hard-coded variables   

    PVC.psfFWHM= 5 % in mm
    PVC.P2Bratio=1.06 % plasma to blood ratio based on 20 volunteer analysis.
    PVC.petMoCo=false; % MR based MoCo only.
    PVC.AC='CT'; % Type of attenuation map - optional
    
    
%% Program start 
    [pvcOutputs] = coreIterativePVC(PVC,r);

