%% This function will perform patlak analysis based on the tissue activity curves and plasma input function (Arterial or image-derived)
% 
% 
%  Author: Lalith Kumar Shiyam Sundar, M.Sc.
%  Date: August 11th, 2018
%  Quantitative Imaging and Medical Physics, Medical University of Vienna.
%  Incase of errors or troubles, contact: lalith.shiyamsundar@meduniwien.ac.at 
%  Usage: [1] patlakInputs.pathOfNiftiPET = path of the nifti PET files. 
%         [2] patlakInputs.pathOfTransformedMRnavigators = path which contains the PET 
%         transformed MR navigators.      
%         [3] patlakInputs.pathOfMRcorrespondenceFile= path of the MR
%         correspondence file, which connects the MR navigators to PET
%         based on least temporal difference.
%         [4] patlakInputs.plasmaGlucose= value for the cold plasma glucose
%         [5] patlakInputs.lumpedConstant= value of the lumped constant,
%         default value is 0.65; 
%         [6] patlakInputs.pathToStore = where to store patlak analysis
%         (Metabolic rate of glucose images will be stored here)
%         [7] patlakInputs.plasmaInputFunction= plasma input function  
%         [8] patlakInputs.BrainMask= BrainMask (binary)
%-------------------------------------------------------------------------%
%                               Program start
%-------------------------------------------------------------------------%


function []=performPatlak(patlakInputs)

%% Hard-coded variables: 
flags.interp=0;

    % Initialization: Hammersmith region values 

    HammerSmithAtlasRegions={'Hippocampus_r_sum';
    'Hippocampus_l_sum';
    'Amygdala_r_sum';
    'Amygdala_l_sum';
    'Ant_TL_med_r_sum';
    'Ant_TL_med_l_sum';
    'Ant_TL_inf_lat_r_sum';
    'Ant_TL_inf_lat_l_sum';
    'G_paraH_amb_r_sum';
    'G_paraH_amb_l_sum';
    'G_sup_temp_cent_r_sum';
    'G_sup_temp_cent_l_sum';
    'G_tem_midin_r_sum';
    'G_tem_midin_l_sum';
    'G_occtem_la_r_sum';
    'G_occtem_la_l_sum';
    'Cerebellum_r_sum';
    'Cerebellum_l_sum';
    'Brainstem_sum';
    'Insula_l_sum';
    'Insula_r_sum';
    'OL_rest_lat_l_sum';
    'OL_rest_lat_r_sum';
    'G_cing_ant_sup_l_sum';
    'G_cing_ant_sup_r_sum';
    'G_cing_post_l_sum';
    'G_cing_post_r_sum';
    'FL_mid_fr_G_l_sum';
    'FL_mid_fr_G_r_sum';
    'PosteriorTL_l_sum';
    'PosteriorTL_r_sum';
    'PL_rest_l_sum';
    'PL_rest_r_sum';
    'CaudateNucl_l_sum';
    'CaudateNucl_r_sum';
    'NuclAccumb_l_sum';
    'NuclAccumb_r_sum';
    'Putamen_l_sum';
    'Putamen_r_sum';
    'Thalamus_l_sum';
    'Thalamus_r_sum';
    'Pallidum_l_sum';
    'Pallidum_r_sum';
    'Corp_Callosum_sum';
    'FrontalHorn_r_sum';
    'FrontalHorn_l_sum';
    'TemporaHorn_r_sum';
    'TemporaHorn_l_sum';
    'ThirdVentricl_sum';
    'FL_precen_G_l_sum';
    'FL_precen_G_r_sum';
    'FL_strai_G_l_sum';
    'FL_strai_G_r_sum';
    'FL_OFC_AOG_l_sum';
    'FL_OFC_AOG_r_sum';
    'FL_inf_fr_G_l_sum';
    'FL_inf_fr_G_r_sum';
    'FL_sup_fr_G_l_sum';
    'FL_sup_fr_G_r_sum';
    'PL_postce_G_l_sum';
    'PL_postce_G_r_sum';
    'PL_sup_pa_G_l_sum';
    'PL_sup_pa_G_r_sum'
    'OL_ling_G_l_sum';
    'OL_ling_G_r_sum';
    'OL_cuneus_l_sum';
    'OL_cuneus_r_sum';
    'FL_OFC_MOG_l_sum';
    'FL_OFC_MOG_r_sum';
    'FL_OFC_LOG_l_sum';
    'FL_OFC_LOG_r_sum';
    'FL_OFC_POG_l_sum';
    'FL_OFC_POG_r_sum';
    'S_nigra_l_sum';
    'S_nigra_r_sum';
    'Subgen_antCing_l_sum';
    'Subgen_antCing_r_sum';
    'Subcall_area_l_sum';
    'Subcall_area_r_sum';
    'Presubgen_antCing_l_sum';
    'Presubgen_antCing_r_sum';
    'G_sup_temp_ant_l_sum';
    'G_sup_temp_ant_r_sum';
    };


if isempty(patlakInputs.lumpedConstant)
    patlakInputs.lumpedConstant=0.65;
end
cd(patlakInputs.pathToStore)
cutOff=patlakInputs.cutOff;
% CMRGlcname='CMRGlcImg_AIF.nii';
% plasmaFileName=dir(['CT_',patlakInputs.subjectID,'*DC*AIF*txt']); % needs to be changed.
CMRGlcname='CMRGlcImg_IDIF_Dixon.nii';
plasmaFileName=dir(['Dixon*',patlakInputs.subjectID,'*native*IDIF*txt']);
if isempty(plasmaFileName)
    plasmaFileName=dir(['PCT*',patlakInputs.subjectID,'*native*IDIF*txt']);
end

%% Create patlak analysis folder 

cd(patlakInputs.pathToStore)
mkdir('Patlak analysis');
pathToPatlakAnalysis=[patlakInputs.pathToStore,filesep,'Patlak analysis'];
cd(pathToPatlakAnalysis);
mkdir('MoCo PET Dixon');
pathToMoCoPET=[pathToPatlakAnalysis,filesep,'MoCo PET Dixon'];

%% Go the folder containing the mr correspondence file 

cd(patlakInputs.pathOfMRcorrespondenceFile);
mrCorrespondenceFile=dir('MR*corres*');
PETtoMRnavLink=readtable(mrCorrespondenceFile.name);
PETtoMRnavLink.Var2=(1:length(PETtoMRnavLink.Var1))';

%% Go the folder containing PET images.
% 
cd(patlakInputs.pathOfNiftiPET);
niftiPETfiles=dir('*nii');
for lp=1:length(niftiPETfiles)
    PETfilesToSort{lp,:}=[patlakInputs.pathOfNiftiPET,filesep,niftiPETfiles(lp).name];
    PETfilesNames{lp,:}=niftiPETfiles(lp).name;
end
sortedPETfiles=natsort(PETfilesToSort);
sortedPETfileNames=natsort(PETfilesNames);

%% Prepare file paths of the PET for rotation 

numberOfUniqueNavigators=unique(PETtoMRnavLink.Var1);
for lp=1:length(numberOfUniqueNavigators)
    indiceOfInterest=PETtoMRnavLink.Var1==numberOfUniqueNavigators(lp);
    mrLinkedPET{lp}=sortedPETfiles(indiceOfInterest);
end

% Go the folder containing the transformed MR navigators in PET space.

cd(patlakInputs.pathOfTransformedMRnavigators)
transformedNiftiMRfiles=dir('Transform*MRnav*nii');
for lp=1:length(transformedNiftiMRfiles)
    mrFilesToSort{lp,:}=transformedNiftiMRfiles(lp).name;
end
sortedMRnavigators=natsort(mrFilesToSort);


% Coregistration of the PET images to a common reference line
% 
CoregInputs.RefImgPath=[patlakInputs.pathOfTransformedMRnavigators,filesep,sortedMRnavigators{1}]; % reference MR navigator
CoregInputs.Interp=1; %Trilinear interpolation
for lp=2:length(numberOfUniqueNavigators)
    CoregInputs.Prefix='MoCo_';
    CoregInputs.SourceImgPath=[patlakInputs.pathOfTransformedMRnavigators,filesep,sortedMRnavigators{numberOfUniqueNavigators(lp)}];
    CoregInputs.MaskImgPath=mrLinkedPET{lp};
%     Coregistration_job(CoregInputs);
    disp(['Performing coregistration based on motion vectors from MR navigator ', num2str(numberOfUniqueNavigators(lp)),'!']);
end

%% storing the motion correction pet in a different folder
% 
% cd(patlakInputs.pathOfNiftiPET)
% movefile('MoCo*',pathToMoCoPET);
% earlyPETframesinAlignment=sortedPETfileNames(PETtoMRnavLink.Var1==numberOfUniqueNavigators(1));
% for lp=1:length(earlyPETframesinAlignment)
%     copyfile(earlyPETframesinAlignment{lp},pathToMoCoPET);
%     cd(pathToMoCoPET)
%     movefile(earlyPETframesinAlignment{lp},['MoCo_',earlyPETframesinAlignment{lp}]);
%     cd(patlakInputs.pathOfNiftiPET);
% end

%Reading the files of the MoCo PET images

cd(pathToMoCoPET)
mocoPETfiles=dir('MoCo*');
for lp=1:length(mocoPETfiles)
    mocoPETtoSort{lp,:}=mocoPETfiles(lp).name;
end
sortedMoCoPET=natsort(mocoPETtoSort);


% Load the brain mask for patlak analysis

cd(patlakInputs.BrainMask)
zeroFile=dir('c1C*nii');
brainMask=spm_read_vols(spm_vol(zeroFile.name));
brainHdr=spm_vol(zeroFile.name);
brainMask=zeros(size(brainMask));
for lp=1:3
    foi=dir(['c',num2str(lp),'C*.nii']);
    brainMask=brainMask+(spm_read_vols(spm_vol(foi.name))>0);
end
brainMask=brainMask>0;

for lp=1:size(brainMask,3);
    brainMask(:,:,lp)=imfill(brainMask(:,:,lp),'holes');
end
brainHdr.fname='Brain-mask.nii';
spm_write_vol(brainHdr,brainMask);
coregT1MRfile=dir('Coreg_T1*');
spm_reslice({coregT1MRfile.name,brainHdr.fname});
brainMask=spm_read_vols(spm_vol('rBrain-mask.nii'));

[PETframingMidPoint]=getPETmidTime();
PETframingMidPoint=round(PETframingMidPoint);

%% Regional patlak only.

if patlakInputs.voxelWise==0
    cd(patlakInputs.plasmaInputFunction); % go to the path containing the plasma input function. 
    plasmaIDIF=readtable(plasmaFileName.name);
    plasmaInputFunction=plasmaIDIF.value_kBq_cc_;
    %plasmaInputFunction=plasmaInputFunction.*1000;
    plasmaInputFunction(plasmaInputFunction<0)=0;
    plasmaInputFunction(isinf(plasmaInputFunction))=0;
    plasmaInputFunction(isnan(plasmaInputFunction))=0;
    %plasmaInputFunctionNonInterp=plasmaInputFunction(PETframingMidPoint);
    integralPlasmaInputFunction=cumsum(plasmaInputFunction); % since the temporal resolution is 1, we can just use cumulative sum (which is actually the integral)
    integralPlasmaInputFunction(plasmaInputFunction==0)=0;
    normIntegralPlasmaInputFunction=integralPlasmaInputFunction./max(integralPlasmaInputFunction);
    integralPlasmaInputFunction(normIntegralPlasmaInputFunction<0.002)=0;
    %xaxisPatlak=integralPlasmaInputFunction./plasmaInputFunctionNonInterp; % normalized time, refer Turku cmrglc patlak webpage.
    xaxisPatlak=integralPlasmaInputFunction./plasmaInputFunction;
    xaxisPatlak=[ones(size(xaxisPatlak)) xaxisPatlak]; % padding it with ones for linear regression.
    xaxisPatlak(isnan(xaxisPatlak))=0;
    xaxisPatlak(isinf(xaxisPatlak))=0;
    xaxisPatlak(xaxisPatlak<0)=0;
    % paths needed for normalization
    cd(patlakInputs.BrainMask)
    coregT1MRfile=dir('Coreg_T1_*.nii');
    pathOfCoregT1MR=[patlakInputs.BrainMask,filesep,coregT1MRfile.name];
    cd(patlakInputs.BrainMask)
    HSatlasFile=dir('Hammers*img')
    pathOfAtlas=[patlakInputs.BrainMask,filesep,HSatlasFile.name];
    deformationFile=dir('iy*nii');
    pathOfDeformationFields=[patlakInputs.BrainMask,filesep,deformationFile.name];
    
    % normalize hammersmith from mni space to native space.
    NHSinputs.DeformationFields=pathOfDeformationFields;
    NHSinputs.HammersmithAtlas=pathOfAtlas;
    normalizeHammerSmith(NHSinputs);
    disp(['Hammersmith atlas normalized to patient space using ',NHSinputs.DeformationFields]);
    normHSfile=dir('Norm*Hammer*img');
    pathToNormHSfile=[NHSinputs.HammersmithAtlas,filesep,normHSfile.name];
    
  
    % Get tissue activity curves from the dynamic PET using the hammersmith atlas.

    % Read in the normalized hammersmith atlas 

    cd(patlakInputs.BrainMask);
    normalizedAtlasFile=dir('Norm*hdr');
    normalizedAtlas=spm_read_vols(spm_vol(normalizedAtlasFile.name));
    normalizedAtlasHdr=spm_vol(normalizedAtlasFile.name);
    voxDim=[abs(normalizedAtlasHdr.private.mat(1,1)) abs(normalizedAtlasHdr.private.mat(2,2)) abs(normalizedAtlasHdr.private.mat(3,3))]

    % Read in the dynamic PET frames
    cd(pathToMoCoPET);
    parfor lp=1:length(sortedMoCoPET)
        dynamicPET{lp}=spm_read_vols(spm_vol(sortedMoCoPET{lp}));
    end
    % calculate regional tissue activity curve
    for olp=1:length(HammerSmithAtlasRegions)
        for ilp=1:length(sortedMoCoPET)
            PETframeOfInterest=dynamicPET{ilp};
            atlasRegionOfInterest=normalizedAtlas==olp;
            regionalTACs{olp}(ilp)=mean(PETframeOfInterest(atlasRegionOfInterest)); %Bq/mL
        end
        TAC=regionalTACs{olp};
        interpTAC=interp1(PETframingMidPoint,TAC,1:3600,'pchip');
        interpTAC=interpTAC';
        yaxisPatlak=interpTAC./plasmaInputFunction;
        %yaxisPatlak=regionalTACs{olp}./plasmaInputFunctionNonInterp';
        yaxisPatlak(isnan(yaxisPatlak))=0;
        yaxisPatlak(isinf(yaxisPatlak))=0;
        yaxisPatlak(yaxisPatlak<0)=0;
        %yaxisPatlak=yaxisPatlak';
        pCoeff=xaxisPatlak\yaxisPatlak;
        plot(xaxisPatlak(:,2),yaxisPatlak,'bo','LineWidth',2,'MarkerFaceColor','b');
        hold on
        plot(xaxisPatlak(:,2),(xaxisPatlak*pCoeff),'r-','LineWidth',2);
        title(HammerSmithAtlasRegions{olp},'FontSize',20);
        ylabel('CT(t)/CP(t)','FontSize',20)
        xlabel('Integral of CP(T)/CP(T)','FontSize',20);
        lgd=legend('Actual data','Fitted data','FontSize', 20);
        lgd.FontSize=20;
        %saveas(gcf,[HammerSmithAtlasRegions{olp},'.png']);
        drawnow
        hold off
        
    end
    
return       
end

%%

% Loading the motion compensation PET images to workspace.
cd(pathToMoCoPET)
parfor lp=1:length(sortedMoCoPET)
    dynamicPET(:,:,:,lp)=(spm_read_vols(spm_vol(sortedMoCoPET{lp,:}))).*brainMask;
    disp(['Loading ',sortedMoCoPET{lp},' image !']);
end

%% Have a dummy nifti header. 

dummyHdr=spm_vol(sortedMoCoPET{end,:});
dummyHdr.fname=CMRGlcname;
dummyHdr.pinfo=[1 0 0]';

%%
% Interpolate plasma input function to a temporal resolution of 1 s.
cd(patlakInputs.plasmaInputFunction); % go to the path containing the plasma input function. 
plasmaIDIF=readtable(plasmaFileName.name);
plasmaInputFunction=plasmaIDIF.value_kBq_cc_;
plasmaInputFunction=plasmaInputFunction.*1000;
%plasmaInputFunctionNonInterp=plasmaInputFunction(PETframingMidPoint)
plasmaInputFunction(isnan(plasmaInputFunction))=0;
plasmaInputFunction(isinf(plasmaInputFunction))=0;
plasmaInputFunction(plasmaInputFunction<0)=0;
integralPlasmaInputFunction=cumsum(plasmaInputFunction); % since the temporal resolution is 1, we can just use cumulative sum (which is actually the integral)
integralPlasmaInputFunction(plasmaInputFunction==0)=0;
normIntegralPlasmaInputFunction=integralPlasmaInputFunction./max(integralPlasmaInputFunction);
integralPlasmaInputFunction(normIntegralPlasmaInputFunction<0.002)=0;
xaxisPatlak=integralPlasmaInputFunction./plasmaInputFunction; % normalized time, refer Turku cmrglc patlak webpage.
indicesToStart=xaxisPatlak>cutOff;xaxisPatlak=xaxisPatlak(indicesToStart);
xaxisPatlak(isnan(xaxisPatlak))=0;% indicesToConsider=xaxisPatlak>cutOff;
xaxisPatlak(isinf(xaxisPatlak))=0;
xaxisPatlak(xaxisPatlak<0)=0;
xaxisPatlak=[ones(size(xaxisPatlak)) xaxisPatlak]; % padding it with ones for linear regression.

% calculate tissue activity curve
[r,c,p]=ind2sub(size(brainMask),find(brainMask>0));
Ki=zeros(size(dynamicPET(:,:,:,1)));V=Ki;
wb=timebar('Performing voxel-wise patlak analysis...'); 
for lp=1:length(r)
    TAC=squeeze(dynamicPET(r(lp),c(lp),p(lp),:));
    TAC(TAC<0)=0;
    if sum(TAC)==0
        continue
    else
        interpTAC=interp1(PETframingMidPoint,TAC,1:3600,'pchip');
        interpTAC=interpTAC';
        yaxisPatlak=interpTAC./plasmaInputFunction;
        yaxisPatlak=yaxisPatlak(indicesToStart);
        yaxisPatlak(isnan(yaxisPatlak))=0;
        yaxisPatlak(isinf(yaxisPatlak))=0;
        yaxisPatlak(yaxisPatlak<0)=0;
        pCoeff = xaxisPatlak\yaxisPatlak ; % solve the system of linear equations.
        pCoeff(isnan(pCoeff))=0;
        pCoeff(isinf(pCoeff))=0;
        pCoeff(pCoeff<0)=0;
        V(r(lp),c(lp),p(lp))=pCoeff(1); %slope
        Ki(r(lp),c(lp),p(lp))=pCoeff(2); %intercept
        timebar(wb,lp/length(r));
%         plot(xaxisPatlak(:,2),yaxisPatlak,'b+','LineWidth',2,'MarkerFaceColor','b ');
%         hold on
%         plot(xaxisPatlak(:,2),(xaxisPatlak*pCoeff),'r-','LineWidth',2);
%         ylabel('CT(t)/CP(t)','FontSize',20)
%         xlabel('Integral of CP(T)/CP(T)','FontSize',20);
%         lgd=legend('Actual data','Fitted data','FontSize', 20);
%         lgd.FontSize=20;
%         pause(0.2)
%         hold off
    end
end
close(wb);
multiplicativeFactor=(patlakInputs.plasmaGlucose./patlakInputs.lumpedConstant);
CMRGlcImg=Ki.*multiplicativeFactor.*6000; % units: (�mol glucose)*(100 g tissue)-1*min-1
figure,imshow3Dfull(CMRGlcImg);
cd(pathToPatlakAnalysis)
delete(patlakInputs.file2Delete)
spm_write_vol(dummyHdr,CMRGlcImg);

% 

