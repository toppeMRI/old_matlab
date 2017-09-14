function [fieldmap,m] = loadfmap(pfile,deltaTE)
% Create fieldmap and magnitude image from Pfile acquired with ~jfnielse/projects/rfovfmri/b0map/scan.tgz.
% See ~jfnielse/projects/rfovfmri/README
% 
% Inputs:
%    pfile
%    deltaTE:  sec
%
% $Id: loadfmap.m,v 1.8 2017/08/08 15:21:03 jfnielse Exp $

if ~exist('deltaTE','var')
	deltaTE = 2.3e-3;    % sec
end

imte1 = recon3dft(pfile,1);                      % [nx ny nz ncoils]
imte2 = recon3dft(pfile,2);
imsos1 = sqrt(sum(abs(imte1).^2,4));
imsos2 = sqrt(sum(abs(imte2).^2,4));
m = sqrt(imsos1.^2+imsos2.^2);

addpath /net/brooks/export/home/jfnielse/matlab/util                  % phasecontrast3d
pc = phasecontrast3d(imte2,imte1);                      % [nx ny nz]
fieldmap = pc / (2*pi*deltaTE);    % Hz

im(fieldmap);

save fieldmap.mat fieldmap m imte1 imte2
