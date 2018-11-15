% zero fill 3D images, Fermi-filtered (in-plane) to reduce ringing
%
% function imout = zeropad(imin, res, fltwidth, type)
%
% imin        [nx ny nz]
% res         [nxtarget nytarget nztarget]      
% fltwidth    Fermi filter transition width (pixels). Default: 2.
% type        'circ' (default) or 'rect'
%

function imout = zeropad(imin, res, fltwidth, type)

if ~exist('type', 'var')
	type = 'circ';
end

import toppe.utils.*

[nx ny nz] = size(imin);
if nx ~= ny
	error('image must be square in x/y');
end

% zero-pad
d = ft3(imin);
dout = zeros(res);
rangex = (res(1)/2+1-nx/2):(res(1)/2+1+nx/2-1);
rangez = (res(3)/2+1-nz/2):(res(3)/2+1+nz/2-1);
dout(rangex, rangex, rangez) = d;

% filter
if ~exist('fltwidth', 'var')
	fltwidth = 2;    % filter transition width (pixels)
end
[fltCirc fltRect] = fermi2d(res(1), size(d,1), fltwidth);
if strcmp('type', 'circ')
	flt = fltCirc;
else
	flt = fltRect;
end
for iz = 1:size(dout,3)
	dout(:,:,iz) = bsxfun(@times, dout(:,:,iz), flt);
end

imout = ift3(dout);

return

