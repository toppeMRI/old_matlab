function im = ift3(D,do3dfft)
%
%	function im = ift3(dat)
%
%	Function does a centered inverse 3DFT of a 3D data matrix.
% 
% $Id: ift3.m,v 1.1 2017/08/02 21:14:20 jfnielse Exp $

if ~exist('do3dfft','var')
	do3dfft = true;
end

if do3dfft
	im = fftshift(ifftn(fftshift(D)));
else
	% don't do fft in 3rd dimension
	for k = 1:size(D,3)
		im(:,:,k) = fftshift(ifftn(fftshift(D(:,:,k))));
	end
end

return;

% EOF
