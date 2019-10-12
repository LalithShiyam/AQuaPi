function [ ROI,C1SliceEnd ] = ROIextractor( ICAsegment,PetrousSlice,CavernousSlice )
%function [ ROI,C1SliceEnd ] = ROIextractor( ICAsegment,PetrousSlice )
% This functions extracts the region of interest "ROI" from the ICA mask,
% which is the internal carotid artery geometry.

if CavernousSlice-PetrousSlice<5
    NoCavernous=1
else
    NoCavernous=0
end

for lp=round(size(ICAsegment,3).*0.50):size(ICAsegment,3)
    MeanIntensity(lp)=mean2(ICAsegment(:,:,lp));
    [~,NumberOfObjects(lp)]=bwlabel(ICAsegment(:,:,lp));
    NormalisedMeanIntensity(lp)=MeanIntensity(lp)./NumberOfObjects(lp); % Mean intensity normalised to the number of objects in each slice.
end

MaybePetrous=PetrousSlice;

StartSlice=MaybePetrous-4; % region around the hill and push the real petrous value
EndSlice=MaybePetrous+4;

PetrousSlice=find(NormalisedMeanIntensity== max(NormalisedMeanIntensity(StartSlice:EndSlice)));


[Peaks,XAxis]=findpeaks(NormalisedMeanIntensity);
CurveForGradCalc=(NormalisedMeanIntensity(1:PetrousSlice));
SlopeValues=diff(CurveForGradCalc);
GradChange=sign(SlopeValues);
GradChange(isnan(GradChange))=0;

for lp=size(GradChange,2):-1:1
    if (GradChange(lp)*-1)>0 | (GradChange(lp)*-1)==0% change made on 20.10.2016
       StartPointVolume=  lp;   
       break
    else
    end
end


%% Extracting the slice below the Cavernous part. Removing the extraction of cavernous - choosing 4 slices after the petrous segment

SecCurveForGradCalc=zeros(size(NormalisedMeanIntensity));
SecCurveForGradCalc(PetrousSlice:CavernousSlice)=NormalisedMeanIntensity(PetrousSlice:CavernousSlice);
SlopeValues=diff(SecCurveForGradCalc);
GradChange=sign(SlopeValues);
GradChange(isnan(GradChange))=0;
InvGradChange=GradChange.*-1

for lp=CavernousSlice-1:-1:PetrousSlice
    if InvGradChange(lp)>0 
        EndPointVolume=lp;
        break
    else 
    end
end

for lp=size(ICAsegment,3):-1:1
    if nnz(ICAsegment(:,:,lp))>0
        temp=lp;
        break
    else
    end
end

if temp==CavernousSlice && NoCavernous
    EndPointVolume=CavernousSlice;
else
end
% SecCurveForGradCalc=zeros(size(NormalisedMeanIntensity));
% SecCurveForGradCalc(PetrousSlice:end)=NormalisedMeanIntensity(PetrousSlice:end);
% SlopeValues=diff(SecCurveForGradCalc);
% GradChange=sign(SlopeValues);
% GradChange(isnan(GradChange))=0;
% 
% 
% for lp=PetrousSlice+1:size(ICAsegment,3)
%     if GradChange(lp)>0 | GradChange(lp)==0 
%         EndPointVolume=lp+3;
%         break
%     else 
%     end
% end



ROI=zeros(size(ICAsegment));
ROI(:,:,StartPointVolume-3:EndPointVolume)=ICAsegment(:,:,StartPointVolume-3:EndPointVolume);
C1SliceEnd=StartPointVolume-3;
end

