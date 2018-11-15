function Dout = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
% Apply low-pass 2D Fermi filter (circular or square) to each slice in a 3D image volume
%
% function Dout = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
%
% D                   2D/3D image volume. Can also be kspace, but in that case set Dtype = 'kspace'.
% fltradius           radius of filter, in pixels
% transitionwidth     Fermi filter roll-off (75%-25% width)
% fltgeom             'circ' (default) or 'rect'
% Dtype               'image' (default) or 'kspace'

import toppe.utils.*

if ~exist('fltgeom', 'var')
	fltgeom = 'circ';
end
if ~exist('Dtype', 'var')
	Dtype = 'image';
end

flt = fermi2d(size(D,1), fltradius, transitionwidth, fltgeom);

for iz = 1:size(D,3)
	if strcmp(Dtype, 'image')
		d = cfftn(D(:,:,iz), 'forward');
	else
		d = D(:,:,iz);
	end
	Dout(:,:,iz) = bsxfun(@times, d, flt);
end

if strcmp(Dtype, 'image')
	Dout = cfftn(Dout, 'inverse');
end

return;
