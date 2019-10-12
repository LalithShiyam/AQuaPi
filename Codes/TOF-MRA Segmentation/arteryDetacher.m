function [ ICA ] = arteryDetacher( VolOhneCOW,ModVolOhneCOW )

nhood = strel('cuboid',[1 3 3]);

%% including a patch work - 12-10-2016 : For eroding the bottom part
img3d=VolOhneCOW;
endslice=round(0.55*(size(VolOhneCOW,3))); temp=zeros(size(VolOhneCOW)); % <-- changed from 0.50 from 0.55 due to vienna patients
temp(:,:,1:endslice)=VolOhneCOW(:,:,1:endslice);toptemp=zeros(size(VolOhneCOW));
toptemp(:,:,endslice+1:end)=VolOhneCOW(:,:,endslice+1:end);
[v1,v2]=VolumePartitioner(toptemp);
v=v1+v2;

ecube=imerode(temp,nhood);
ecube(:,:,endslice-5:end)=img3d(:,:,endslice-5:end);% 21-10-2016
[erodedICA]=ECAprune(ecube); % <-- this is creating a lot of issues. use a better algorithm, icaextractor replaced by ecaprune

%% Dilate the internal carotid artery

[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 3; 
[edilate]=imdilate(erodedICA(:,:,1:endslice),nhood); 

tempvol=zeros(size(VolOhneCOW));
tempvol(:,:,1:endslice)=edilate(:,:,1:endslice);
tempvol(:,:,endslice+1:end)=v(:,:,endslice+1:end);

%% Get only the ICA

[ICA]=tempvol.*ModVolOhneCOW;
ICA=ICA>0;[ICA1,ICA2]=VolumePartitioner(ICA);ICA=ICA1+ICA2;
end

