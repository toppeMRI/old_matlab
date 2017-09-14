function [ims imsos d]= recon3dft(pfile,echo,clim)
% function [ims imsos]= recon3dft(pfile,echo)
%
% Recon 3D spin-warp image
%
% Input:
%  echo:         echo to recon (1 or 2)
%
% Output:
%  ims:          [nx ny nz ncoils]    
%
% $Id: recon3dft.m,v 1.10 2017/08/08 15:21:03 jfnielse Exp $

if ~exist('echo','var')
	echo = 1;
end

% load raw data
addpath /net/brooks/export/home/jfnielse/github/toppe/matlab/lib/GE   % loadpfile.m
d = loadpfile(pfile);   % int16. [ndat ncoils nslices nechoes nviews] = [ndat ncoils nz 2 ny]
d = permute(d,[1 5 3 2 4]);         % [ndat ny nz ncoils nechoes]
d = double(d);

%if(mod(size(d,3),2))
%	d = d(:,:,2:end,:,:);  % throw away dabslice = 0
%end

% get flat portion of readout
addpath /net/brooks/export/home/jfnielse/github/toppe/matlab/lib      % readmod.m
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod('readout.mod');
nramp = 0; %15;  % see mat2mod.m
nbeg = paramsint16(3) + nramp;  
nx = paramsint16(4);  % number of acquired data samples per TR
decimation = paramsint16(10);
d = d(nbeg:(nbeg+nx-1),:,:,:,:);     % [nx*125/oprbw ny nz ncoils nechoes]

% recon 
for coil = 1:size(d,4)
	fprintf(1,'recon coil %d\n', coil);
	imstmp = ift3(d(:,:,:,coil,echo));
	ims(:,:,:,coil) = imstmp(end/2+((-nx/decimation/2):(nx/decimation/2-1))+1,:,:);
end

% flip x
%ims = flipdim(ims,1);

% display root sum-of-squares image
imsos = sqrt(sum(abs(ims).^2,4)); 
%figure; im('blue0',imsos,[0 1.3]);
if exist('clim','var')
	figure; im(permute(imsos,[2 1 3]),clim);
else
	figure; im(permute(imsos,[2,1,3]));
end

return;
