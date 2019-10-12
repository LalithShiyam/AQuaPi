%% Purpose
% This is the core program which is responsible for the automatic petrous
% and cervical segmentation of the internal carotid arteries. The
% segmentation is based on Gibo et.al
%% References
% Gibo H, Lenkey C and Rhoton AL. Microsurgical anat- omy of the supraclinoid portion of the internal carotid artery. J Neurosurg 1981; 55: 560–574.  

%%
%% Author information
% Lalith Kumar Shiyam Sundar, 
% Quantitative Imaging and Medical Physics, Medical University of Vienna

%% Inputs
% coreSegInputs.pathOfAngio - Physical path of the 3D time-of-flight MR
% angiography DICOM images. 
% This will be provided by the "wrapperForSegmentation.m"

%% Program start

function []=coreSegmentation(coreSegInputs)

%% Hard-coded variables: change it according to your needs.

upperSliceLimit = 175;
quantileVal     = 0.985; % play around with this value for over or under segmentation
wrapPrevention  = 50; 

%% Logic start
%%
% 
% 1. Carotid vasculature segmentation
% 

% Read in the DICOM images

orgVolume=MPRage(coreSegInputs.pathOfAngio);
orgVolume(:,:,upperSliceLimit:end)=zeros; % nulling the values beyond this point.

% Preparing the struct parameter for entire carotid vasculature
% segmentation.

segCVinputs.imgVolume=orgVolume;
segCVinputs.qVal=quantileVal;
[ICA]=segCV(segCVinputs); % Carotid vasculature segmentation
%%
%
% 2. Morphometric analysis

% Inputs
IFAinputs.orgVolume=orgVolume;
IFAinputs.bVol=ICA;
[IFAoutputs] = intensityFeatureAnalysis(IFAinputs); % feature calculation

%%
%
% 3. Circle of willis (COW) detection and removal
CoWDetectorInputs.bVol=cropBorders(ICA,wrapPrevention);
CoWDetectorInputs.curveOfInterest=IFAoutputs.curveOfInterest;
[icaWithoutCOW] = CoWdetector( CoWDetectorInputs);

%%
%
% 4. Major axis length calculation

[ prepICA ] = MALOrientationCalc( icaWithoutCOW ); % calculate major axis length

%%
%
% 6. Remove overlay between ICA and ECA for ICA isolation.
% 
remOLinputs.orgVol=orgVolume;
remOLinputs.prepICA=prepICA;
[modifiedICA]=remOverlay(remOLinputs);% removes the overlay between ICA and ECA

%%
%
% 7. Detaching the ICA and ECA
%
[ICA]=arteryDetacher(modifiedICA,prepICA); 

%%
%
% 8. Getting petrous slices

[ vol1,vol2 ] = VolumePartitioner( ICA ); % separates into two volumes
 [petrousSlice1]=PetrousFinder(vol1,orgVolume); % calculates the petrous slice for caroitd 1
[petrousSlice2]=PetrousFinder(vol2,orgVolume); % calculates th epetrous slice for the carotid 2
[ cavernousSlice1 ] = CavernousSliceDetector( petrousSlice1,vol1 ); % calculates the Cervical segment for the carotid 1
[ cavernousSlice2 ] = CavernousSliceDetector( petrousSlice2,vol2 ); % calculates the Cervical segment for the carotid 2
[ PreROI1,C1SliceEnd1 ] = ROIextractor( vol1,petrousSlice1,cavernousSlice1 ); % extracts the ROI
[ PreROI2,C1SliceEnd2 ] = ROIextractor( vol2,petrousSlice2,cavernousSlice2 );
PreROI=PreROI1+PreROI2;
[ROI1,ROI2]=VolumePartitioner(PreROI);
c2Seg=ROI1+ROI2; % petrous segment
[temp1,temp2]=VolumePartitioner(ICA);

%%
%
% 9. Getting cervical slices

% Extraction of the cervical region.

[ C1SegmentICA1 ] = C1Extractor( vol1,C1SliceEnd1,orgVolume );
[ C1SegmentICA2 ] = C1Extractor( vol2,C1SliceEnd2,orgVolume );
[ FC1SegmentICA1 ] = FinalCleanUp( C1SegmentICA1,vol2);
[ FC1SegmentICA2 ] = FinalCleanUp( C1SegmentICA2,vol2);


C1SegmentICA=FC1SegmentICA1+FC1SegmentICA2;
[C1seg1,C1seg2]=VolumePartitioner(C1SegmentICA);
C1Seg=C1seg1+C1seg2; % Cervical segment



whereToStore=coreSegInputs.path2StoreSeg;
DICOMwriter([coreSegInputs.patientCode,'_','Petrous'],C2Seg,'C2Seg',whereToStore,1);
DICOMwriter([coreSegInputs.patientCode,'_','Cervical'], C1Seg,'C1Seg',whereToStore,1);
figure,VolumeVisualizer(C1Seg,'r');title('Cervical')
figure,VolumeVisualizer(C2Seg,'r');title('Petrous')
end




%% Local programs

%%
% 
% 1. For Segmenting the entire carotid vasculature.
% 
function [ bVol ] = segCV( segCVinputs )

% passing inputs on to local variables.

orgVol=segCVinputs.imgVolume;
qVal=segCVinputs.qVal;

% preliminary segmentation

qMask=(orgVol>quantile(orgVol(:),qVal));

% Selecting a particular image slice where the ICA is bright for automatic
% seeded region growing algorithm

orgVol2D=orgVol(:,:,round(size(orgVol,3)*0.5));
orgVol2D=double(orgVol2D); 

% Thresholding for getting the seed-points

threshold=round(0.75*(max(orgVol2D(:))));
BM=orgVol2D>threshold;
skelImg=bwmorph(BM,'skel',inf);
seedImg=bwmorph(skelImg,'shrink');
[r,c]=find(seedImg==1);

% Region growing algorithm

% Getting the automatically selected seeds as matrix
seedMatrix=[r,c];
seedMatrix=[seedMatrix,round(size(orgVol,3)*0.5)*ones(size(seedMatrix,1),1)];
bVol=zeros(size(orgVol));       

for lp=1:size(seedMatrix,1)
    [~,maskOfInt]=regionGrowing(qMask,seedMatrix(lp,:));
    bVol=bVol+maskOfInt;
end
bVol=(bVol>=1);
[V1,V2]=VolumePartitioner(bVol);
clear bVol
bVol=V1+V2;
bVol=(bVol>=1);
end
    
%%
%
% 2. Morphometric analysis

function [IFAoutputs] = intensityFeatureAnalysis( IFAinputs )

% Passing inputs to local variables.
orgVol=IFAinputs.orgVolume;
bVol=IFAinputs.bVol;
gVol=orgVol.*bVol;
endPt=size(gVol,3);

% Initialisation for optimisation purposes

NumberOfObjects=zeros(1,size(gVol,3));
MeanIntensity=zeros(1,size(gVol,3));
NormalisedMeanIntensity=zeros(1,size(gVol,3));
MInonzero=zeros(1,size(gVol,3));
NMInonzero=zeros(1,size(gVol,3));

% Calculating the Mean intensity and the normalised mean intensity.

for lp=1:endPt
    MeanIntensity(lp)=mean2(gVol(:,:,lp));
    MInonzero(lp)=mean(nonzeros(gVol(:,:,lp)));
    [~,NumberOfObjects(lp)]=bwlabel(bVol(:,:,lp));
    NormalisedMeanIntensity(lp)=MeanIntensity(lp)./NumberOfObjects(lp); % Mean intensity normalised to the number of objects in each slice.
    NMInonzero(lp)=MInonzero(lp)./NumberOfObjects(lp);
end

% Getting only the features for the top part of the volume 

curveOfInterest=NormalisedMeanIntensity;
IFAoutputs.curveOfInterest=curveOfInterest;
end

%%
%
% 3. Cropping the volume to escape wrap-around artifacts.

%% cropping the image volume.
function [croppedVolume]=cropBorders(bVol,wrapPrevention)
cropFromBorder=wrapPrevention; % pixels
croppedVolume=zeros(size(bVol));
mask=zeros(size(bVol));
for lp=1:size(bVol,3)
    maskZero=ones(size(bVol(:,:,lp)));
     maskZero(1:cropFromBorder,:)=0;
     maskZero(:,1:cropFromBorder)=0;
     maskZero(:,(size(maskZero,2)-cropFromBorder):end)=0;
     maskZero((size(maskZero,1)-cropFromBorder):end,:)=0;
     mask(:,:,lp)=maskZero;
     croppedVolume(:,:,lp)=mask(:,:,lp).*bVol(:,:,lp);
end
end

%%
%
% 4. Removing the circle of willis
% 
function [icaWithoutCOW] = CoWdetector( CoWDetectorInputs)
curveOfInterest=CoWDetectorInputs.curveOfInterest;
bVol=CoWDetectorInputs.bVol;
limit=round(0.65*length(curveOfInterest));
temp=zeros(size(curveOfInterest));
tempNMI=zeros(size(curveOfInterest));
tempNMI(limit:end)=curveOfInterest(limit:end);
smoothTempNMI=sgolayfilt(tempNMI,4,15);

%% Major axis length - might be the salvation - check.
[~,endInfo]=StartEndImgInfo(bVol);
for lp=limit:endInfo
    G=bwlabel(bVol(:,:,lp));
    FillArea=regionprops(G,'FilledArea');
    threshVF=double(round(max([FillArea.FilledArea])/5));   
    if  isempty(threshVF)
        threshVF=1
    end
    
    BW2 = bwareaopen(bVol(:,:,lp),threshVF);
    if nnz(BW2)==0
        temp(lp)=0;
    else
        [m1,~]=VolumePartitioner(BW2);
        Stats=regionprops(m1,'MajorAxisLength','MinorAxisLength');
        temp(lp)=Stats.MajorAxisLength;
    end
    
end

% Major/Minor axis quantification for every object in each slice.

for lp=1:endInfo
    G=bwlabel(bVol(:,:,lp));
    FillArea=regionprops(G,'FilledArea');
    threshVF=double(round(max([FillArea.FilledArea])/5));   
    if  isempty(threshVF)
        threshVF=1
    end
    
    BW2 = bwareaopen(bVol(:,:,lp),threshVF);
    if nnz(BW2)==0
        majorMinor(lp)=0;
    else
        [m1,m2]=VolumePartitioner(BW2);
        cc = bwconncomp(BW2);
        NumberOfObjects=cc.NumObjects;
        Stats=regionprops(m1,'MajorAxisLength','MinorAxisLength');
        majorMinor(lp)=Stats.MajorAxisLength/Stats.MinorAxisLength;
    end
end

smoothMajorMinor=sgolayfilt(majorMinor,4,15);
smoothMajorMinor(endInfo+1:size(bVol,3))=0;

% Morphometric analysis

curveOfInterest=temp.*tempNMI;
smoothCOI=sgolayfilt(curveOfInterest,4,7); % filter length changed from 15 to 7
smoothCOI(isnan(smoothCOI))=0;
curveOfInterest=smoothCOI.*smoothMajorMinor;
[Peaks,XAxis]=findpeaks(curveOfInterest);
[~,SortedIdx]=sort(Peaks,'descend');
sortedXAxis=XAxis(SortedIdx);
sortedXAxis=sortedXAxis(1:4);

for lp=1:length(sortedXAxis)
    conncomp=bwconncomp(bVol(:,:,sortedXAxis(lp)));
    [~,idx]=max(cellfun('size',conncomp.PixelIdxList,1));
    BM=zeros(size(bVol(:,:,sortedXAxis(lp))));
    BM(conncomp.PixelIdxList{idx})=1;
    Stats=regionprops(BM,'MajorAxisLength','Orientation');
    MALMatrix(lp,:)=[sortedXAxis(lp),Stats.MajorAxisLength];
    OrientationMatrix(lp,:)=[sortedXAxis(lp),Stats.Orientation];
end

% geometry box --> remove the anomoly of COW.

for lp=1:length(MALMatrix)
    [t1,t2]=VolumePartitioner(bVol(:,:,MALMatrix(lp,1)));
    stats=regionprops(t1,'BoundingBox')
    widthOfRectangle(lp,1)=MALMatrix(lp,1)
    widthOfRectangle(lp,2)=stats.BoundingBox(4) % we are getting the width and 3 corresponds to the width section.
end

for lp=1:length(widthOfRectangle)
    bw=bVol(:,:,widthOfRectangle(lp,1));
    [t1,t2]=VolumePartitioner(bw);stats=regionprops(t1,'BoundingBox')
    figure,imshow(t1)
    hrectangle=imrect(gca,[stats.BoundingBox]);
    maskimg=hrectangle.createMask();OrientationOfPetrous(lp,1)=widthOfRectangle(lp,1);
    temp=regionprops(maskimg,'Orientation');
    OrientationOfPetrous(lp,2)=temp.Orientation;
end

idx=(OrientationOfPetrous(:,2)==90); 
probCavernous=OrientationOfPetrous(idx,1);
probCavernous=sort(probCavernous,'ascend');
probPetrousRegion=OrientationOfPetrous(~idx,1);
probPetrousRegion=sort(probPetrousRegion,'ascend');

% calculating the slopes 

curveForGradCalc=zeros(size(CurveOfInterest));
curveForGradCalc(probCavernous:end)=CurveOfInterest(probCavernous:end);
slopeValues=diff(curveForGradCalc);
gradChange=sign(slopeValues);
gradChange(isnan(gradChange))=0;

for lp=probCavernous:size(gradChange,2)
    if GradChange(lp)< 0
       endPointVolume= lp;
    else
        break
    end
end

% Getting the ICA without the Circle of willis

%% Extracting the two largest connected components 

icaWithoutCOW=zeros(size(bVol));
tempVol=icaWithoutCOW;
tempVol=bVol(:,:,1:(endPointVolume+1));
conncomp=bwconncomp(tempVol);
[~,i4]=max(cellfun('size',conncomp.PixelIdxList,1));
icaWithoutCOW(conncomp.PixelIdxList{i4})=1;
conncomp.PixelIdxList{i4}=[];
[~,i4]=max(cellfun('size',conncomp.PixelIdxList,1));
icaWithoutCOW(conncomp.PixelIdxList{i4})=1;

end
%%
%
% 5. Major axis length calculation

function [ prepICA ] = MALOrientationCalc( bVol )

[startinfo,endinfo]=StartEndImgInfo(bVol);

for lp=startinfo:size(bVol,3)
    if nnz(bVol(:,:,lp))==0
        EndSlice=lp;
        break
    else
        if lp==size(bVol,3) 
            EndSlice=size(bVol,3);
        else
        end
    end
end

StartSlice=round(0.80*EndSlice);

for lp=StartSlice:EndSlice
    if nnz(bVol(:,:,lp))==0
       ObjOrientation(lp)=0;
       ObjMajorAxisLength(lp)=0;
    else
        conncomp=bwconncomp(bVol(:,:,lp));
        [~,idx]=max(cellfun('size',conncomp.PixelIdxList,1));
        BM=zeros(size(bVol(:,:,lp)));
        BM(conncomp.PixelIdxList{idx})=1;
        Stats=regionprops(BM,'MajorAxisLength','Orientation');
        ObjOrientation(lp)=Stats.Orientation;
        ObjMajorAxisLength(lp)=Stats.MajorAxisLength;
    end
end


%% Finding Peaks for the Major axis length vector.

[peaks,locs]=findpeaks(ObjMajorAxisLength);
[Val,Idx]=sort(peaks,'descend');
MALSOI=locs(Idx);MALSOI=MALSOI(1);

%% Finding the orientation of the object in the slice of interest.

if ObjOrientation(MALSOI) < 15 && ObjOrientation(MALSOI) > -15
    SliceOfInterest= MALSOI;
    ModVol=zeros(size(bVol));
    ModVol(:,:,MALSOI)=zeros(size(bVol(:,:,MALSOI)));
    ModVol(:,:,1:MALSOI-1)=bVol(:,:,1:MALSOI-1);
    [fv1,fv2]=VolumePartitioner(ModVol);
    prepICA=fv1+fv2;
else
    prepICA=bVol;
end

end
%%
%
% 6. Remove overlay between ICA and ECA for ICA isolation.

function [modifiedICA]=remOverlay(remOLinputs)
orgVolume=remOLinputs.orgVolume;
prepVol=remOLinputs.prepICA;
temp=orgVolume.*prepVol;

endSlice=round(0.40*size(orgVolume,3));

for lp=1:endSlice
    I_eq=temp(:,:,lp);  
    bw = imbinarize(I_eq, graythresh(I_eq));
    bw2 = imfill(bw,'holes');
    bw3 = imopen(bw2, ones(5,5));
    bw4 = bwareaopen(bw3, 30);
    bw4_perim = bwperim(bw4);
    overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
    mask_em = imextendedmax(I_eq, 30); %proposed was 30 and changed to 40 and now its back to 30 again.
    mask_em = imclose(mask_em, ones(5,5));
    mask_em = imfill(mask_em, 'holes');
    G=bwlabel(mask_em);
    FillArea=regionprops(G,'FilledArea');
    threshVF=double(round(max([FillArea.FilledArea])/5));   
    if  isempty(threshVF)
        threshVF=1
    end
    mask_em = bwareaopen(mask_em, threshVF);
    overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
    I_eq_c = imcomplement(I_eq);
    I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
    L(:,:,lp) = watershed(I_mod); 
end
modifiedICA=zeros(size(temp));
modifiedICA(:,:,1:endSlice)=prepVol(:,:,1:endSlice);
modifiedICA(L == 0) = 0;
modifiedICA(:,:,endSlice:end)=prepVol(:,:,endSlice:end);

end

