function [ CavernousSlice ] = CavernousSliceDetector( PetrousSlice,ICAsegment )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for lp=PetrousSlice+7:size(ICAsegment,3) %bad patch added 18-10-2016 , to cover the hill part. (petrousslice+3)
    MI(lp)=mean2(ICAsegment(:,:,lp));
    [~,NumberOfObjects(lp)]=bwlabel(ICAsegment(:,:,lp));
    NMI(lp)=MI(lp)./NumberOfObjects(lp); % Mean intensity normalised to the number of objects in each slice.
    NMI(isnan(NMI))=0;
end

for lp=size(ICAsegment,3):-1:1
    if nnz(ICAsegment(:,:,lp)>0)
       EndVal=lp;
       break
    else
    end
end

[Pks,Locs]=findpeaks(NMI); 
if isempty(Locs)
    CavernousSlice=EndVal;
else
    if length(Locs)==1
       CavernousSlice=Locs;
    else
        [Vals,Idx]=sort(Pks,'Descend')
        SortLocs=Locs(Idx);
        %% write a code snippet which handles a condition where there is only two minimas

        %% need an efficient code. : patchwork applied at 15.06.2016

        if length(SortLocs)==2
        else
            SortLocs=SortLocs(1:3);
        end

        for lp=1:length(SortLocs)
            temparea=regionprops(ICAsegment(:,:,SortLocs(lp)),'Area')
            AreaMatrix(lp,:)=[SortLocs(lp),temparea.Area];
            temporientation=regionprops(ICAsegment(:,:,SortLocs(lp)),'Orientation')
            OrientationMatrix(lp,:)=[SortLocs(lp),temporientation.Orientation];
            tempMal=regionprops(ICAsegment(:,:,SortLocs(lp)),'MajorAxisLength')
            MALMatrix(lp,:)=[SortLocs(lp),tempMal.MajorAxisLength];
        end

        % LogIdx=(MALMatrix(:,1)==PetrousSlice);
        % 
        % MALpetrous=MALMatrix(LogIdx,2);
        % 
        % CavMALMatrix=MALMatrix(~LogIdx,:);
        % 
        % % CavMALMatrix(CavMALMatrix(:,2)>MALpetrous)=[];
        % % 
        % % DiffVec=(MALpetrous-CavMALMatrix(:,2));
        % % 
        % % IdxCavernous=(DiffVec==min(DiffVec));
        % % 
        % % CavernousSlice=CavMALMatrix(IdxCavernous,1);
        % 
        % idx=find(CavMALMatrix(:,2)==max(CavMALMatrix(:,2)));
        % CavernousSlice=CavMALMatrix(idx,1);
        % 
        % CavIdx=find(MALMatrix(:,2)==max(MALMatrix(:,2)))
        % CavernousSlice=MALMatrix(CavIdx,1);

      %  CavernousSlice=find(NMI==max(NMI));

        for lp=1:size(OrientationMatrix,1)
            if OrientationMatrix(lp,2)>-20 && OrientationMatrix(lp,2)<20
               %OrientationMatrix(lp,1)=0;
               MALMatrix(lp,2)=0;
            else
            end
        end
     %CavernousSlice=find(NMI==max(NMI));
      temp=(MALMatrix(:,2)==max(MALMatrix(:,2)));CavernousSlice=MALMatrix(temp) 
    end
end



end

