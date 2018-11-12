function imout = zeropad(imin, res, fltwidth)
% zeropad 3D images, Fermi-filtered (in-plane) to reduce ringing
%
% imin        [nx ny nz]
% res         [nxtarget nytarget nztarget]      
% fltwidth    Fermi filter transition width (pixels). Default: 2.
%
% $Id: zeropad.m,v 1.4 2018/11/12 16:30:49 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/zeropad.m,v $

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
flt = sub_fermi2d(res(1), size(d,1), fltwidth);
for iz = 1:size(dout,3)
	dout(:,:,iz) = bsxfun(@times, dout(:,:,iz), flt);
end

imout = ift3(dout);

return


function f= sub_fermi2d(sz,radius,width)
%  usage ...  fermi2d(sz,radius,width)
%  sz - matrix size, radius - radius of window, width - transition width
%  negative width will give "standard width" = 20*radius/128
%
% Doug Noll

%
% [x,y]= meshdom(-64:1:63,  63:-1:-64);
% f= 1 ./ (1 + exp((sqrt(x.^2 + y.^2) - radius) ./ (10*steepness/256)));
if width < 0,
   width = 20*radius/128;
end
i = sqrt(-1);
cent = sz/2 + 1;
x= (1:sz);
y= (1:sz)';
X= ones(size(y))*x;
Y= y*ones(size(x));
clear x y;
R =  abs( X-cent + (Y-cent).*i );
clear X Y;
f = 1 ./ (1 + exp( (R - radius) ./ width ));
clear R;
