function [ ICA ] = FinalCleanUp( VolOhneCOW,ModVolOhneCOW)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Creating a 3d Structuring element
nhood = strel('cube',[3] );
%  [xx,yy,zz] = ndgrid(-5:5);  
%   nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.8; % give a variable formula instead of a value, 2 for denmark and 1.4 for vienna <---changed from 2 to 0.9 due to the loops.

%% Erode the Volume without the circle of willis.

%% including a patch work - 12-10-2016 : For eroding the bottom part
img3d=VolOhneCOW;
[~,endslice]=StartEndImgInfo(img3d);
temp=zeros(size(VolOhneCOW)); % <-- changed from 0.50 from 0.55 due to vienna patients
temp(:,:,1:endslice)=VolOhneCOW(:,:,1:endslice);
[v1,v2]=VolumePartitioner(temp);
v=v1+v2;
%---------------------------
%ecube=imerode(VolOhneCOW,nhood);
ecube=imerode(temp,nhood);
%patch--------------
ecube(:,:,endslice-5:end)=img3d(:,:,endslice-5:end);% 21-10-2016
[bigpart,~]=VolumePartitioner(ecube);
    
%% Dilate the internal carotid artery

  [xx,yy,zz] = ndgrid(-5:5);
 nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2.5; % 2.6 for denmark and 2.5 for vienna
%nhood=strel('cube',4);
[edilate]=imdilate(bigpart(:,:,1:endslice),nhood); %patch work

%% patch 12-10-2016

%edilate(:,:,endslice-5:end)=img3d(:,:,endslice-5:end);

% temp patch 12-10-2016

tempvol=zeros(size(VolOhneCOW));
tempvol(:,:,1:endslice)=edilate(:,:,1:endslice);
tempvol(:,:,endslice+1:end)=v(:,:,endslice+1:end);

%% Get only the ICA

%[ICA]=edilate.*VolOhneCOW;
[ICA]=tempvol.*ModVolOhneCOW;
ICA=ICA>0;[ICA1,ICA2]=VolumePartitioner(ICA);ICA=ICA1+ICA2;
end

