function [ PetrousSlice ] = PetrousFinder( ICAsegment,OriginalVolume )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

StartPoint=round(0.50*size(ICAsegment,3));
TempWorkVol=OriginalVolume.*ICAsegment;

for lp=size(ICAsegment,3):-1:1
    if nnz(ICAsegment(:,:,lp)>0)
        EndPoint=lp;
        break
    else
        continue
    end
end

for lp=StartPoint:EndPoint
    conncomp=bwconncomp(ICAsegment(:,:,(lp)));
    [~,idx]=max(cellfun('size',conncomp.PixelIdxList,1));
    BM=zeros(size(ICAsegment(:,:,lp)));
    BM(conncomp.PixelIdxList{idx})=1;
    Stats=regionprops(BM,'MajorAxisLength','Orientation','Area');
    AreaMatrix(lp)=[Stats.Area];
    OrientationMatrix(lp)=[Stats.Orientation];
    MALMatrix(lp)=[Stats.MajorAxisLength];
end

%
MALMatrix=MALMatrix';
[Val,Idx]=findpeaks(MALMatrix);
WorkMatrix=[Val,Idx];
[SortVal,SortIdx]=sort(WorkMatrix(:,1),'descend');
MayBePetrous=Idx(SortIdx);
MayBePetrous=MayBePetrous(1:2);

%% trial : weird patch 20-10-2016 : applied to handle the branch of circle of willis.
% 
% for lp=1:length(MayBePetrous)
%     conncomp=bwconncomp(ICAsegment(:,:,MayBePetrous(lp)))
%     if conncomp.NumObjects>4
%         MayBePetrous(lp)=[]
%     else
%     end
% end

%%  Gradient check along y direction - to see if its a feature 

% for lp=1:length(MayBePetrous)
%     [~,gy{lp}]=imgradientxy(ICAsegment(:,:,MayBePetrous(lp)));
% end
% 
% for lp=1:length(gy)
%     %Stats=regionprops(gy{lp},'Area','MajorAxisLength');
%     %PetrousProp{lp}=[MayBePetrous(lp),Stats.Area,Stats.MajorAxisLength];
%     petsum(lp)=nnz(gy{lp});
% end
% 
% PetrousSlice=MayBePetrous(find(petsum==max(petsum)));


%PetrousSlice=find(NMI==PetrousMean)

%% Geometry box to find petrous 


for lp=1:length(MayBePetrous)
    [t1,t2]=VolumePartitioner(ICAsegment(:,:,MayBePetrous(lp)));
    stats=regionprops(t1,'BoundingBox')
    widthOfRectangle(lp,1)=MayBePetrous(lp,1)
    widthOfRectangle(lp,2)=stats.BoundingBox(4) % we are getting the width and 4 corresponds to the width section.
end
    
% [~,SortedIdx]=sort(widthOfRectangle(:,2),'descend');
% PetrousCandidates=widthOfRectangle(SortedIdx,1);    
% PetrousCandidates=PetrousCandidates(1:2);

for lp=1:length(widthOfRectangle)
    bw=ICAsegment(:,:,widthOfRectangle(lp,1));
    [t1,t2]=VolumePartitioner(bw);stats=regionprops(t1,'BoundingBox')
    figure,imshow(t1)
    hrectangle=imrect(gca,[stats.BoundingBox]);
    maskimg=hrectangle.createMask();OrientationOfPetrous(lp,1)=widthOfRectangle(lp,1);
    temp=regionprops(maskimg,'Orientation');
    OrientationOfPetrous(lp,2)=temp.Orientation;
end

idx=(OrientationOfPetrous(:,2)==90); 
ProbCavernous=OrientationOfPetrous(idx,1);
PetrousSlice=OrientationOfPetrous(~idx,1);


end

