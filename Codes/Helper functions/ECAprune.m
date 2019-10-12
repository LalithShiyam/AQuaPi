function [ ICA] = ECAprune( erodedICA )

[fv1,fv2]=VolumePartitioner(erodedICA);

% cleaning for small pixels <-- comment to self : why to do this, we
% already used a connected component analysis, (...seriously).
[startinfo1,endinfo1]=StartEndImgInfo(fv1);
[startinfo2,endinfo2]=StartEndImgInfo(fv2);

Volfirst=zeros(size(fv1));
Volsecond=zeros(size(fv2)); 
endval1=round(0.6*endinfo1);endval2=round(0.6*endinfo2);

%% for first volume

for lp=1:endval1;
    L=bwlabel(fv1(:,:,lp));
    FillArea=regionprops(L,'FilledArea');
    threshVF=double(round(max([FillArea.FilledArea])/10));
    if  isempty(threshVF)
        threshVF=1
    end
    Volfirst(:,:,lp)=bwareaopen(fv1(:,:,lp),threshVF);% it was ten before, changed to 8 <-- bad idea, make it dynamic.
end

%%for second volume
for lp=1:endval2
    S=bwlabel(fv2(:,:,lp));
    FillArea2=regionprops(S,'FilledArea');
    threshVS=double(round(max([FillArea2.FilledArea])/10)); 
    if isempty(threshVS)
        threshVS=1
    end
    Volsecond(:,:,lp)=bwareaopen(fv2(:,:,lp),threshVS);
end

%

[Volfirst,~]=VolumePartitioner(Volfirst);
[Volsecond,~]=VolumePartitioner(Volsecond);

Volfirst(:,:,endval1:end)=fv1(:,:,endval1:end);
Volsecond(:,:,endval2:end)=fv2(:,:,endval2:end);

[Volfirst,~]=VolumePartitioner(Volfirst);
[Volsecond,~]=VolumePartitioner(Volsecond);


%% Processing block for the first volume - volfirst (one of the carotids)

%finding the branching point

for lp=1:round(0.50*size(Volfirst,3)) % Searching only the first half for bifurcations.
    [~, numberOfObject] = bwlabel(Volfirst(:,:,lp));
    if numberOfObject == 1 && lp==round(0.50*size(Volfirst,3))
        startSlice=1;
        ICA1=Volfirst;
        break
    else
        if numberOfObject>1
            startSlice=lp;
            break
        else
            continue
        end
    end    
end

% finding the thickest part.
ICA1=zeros(size(Volfirst));
tempVolFirst=zeros(size(Volfirst));
tempVolFirst(:,:,startSlice+1:end)=Volfirst(:,:,startSlice+1:end);
conncomp=bwconncomp(tempVolFirst);
[~,i4]=max(cellfun('size',conncomp.PixelIdxList,1));
ICA1(conncomp.PixelIdxList{i4})=1;
ICA1(:,:,1:startSlice+1)=Volfirst(:,:,1:startSlice+1);

%% Processing block for the second volume
for lp=1:round(0.50*size(Volsecond,3)) % Searching only the first half for bifurcations.
    [~, numberOfObject] = bwlabel(Volsecond(:,:,lp));
    if numberOfObject == 1 && lp==round(0.50*size(Volsecond,3))
        startSlice=1;
        ICA2=Volsecond;
        break
    else
        if numberOfObject>1
            startSlice=lp;
            break
        else
            continue
        end
    end    
end

% finding the thickest part.
ICA2=zeros(size(Volsecond));
tempVolsecond=zeros(size(Volsecond));
tempVolsecond(:,:,startSlice+1:end)=Volsecond(:,:,startSlice+1:end);
conncomp=bwconncomp(tempVolsecond);
[~,i4]=max(cellfun('size',conncomp.PixelIdxList,1));
ICA2(conncomp.PixelIdxList{i4})=1;
ICA2(:,:,1:startSlice+1)=Volsecond(:,:,1:startSlice+1);

%% adding both the volumes

ICA=ICA1+ICA2;
end
%figure,isosurface(ICA,0),axis equal

