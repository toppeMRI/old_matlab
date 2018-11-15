function [fCirc, fRect] = fermi2d(sz,radius,width)
%  usage ...  fermi2d(sz,radius,width)
%  sz - matrix size, radius - radius of window, width - transition width
%  negative width will give "standard width" = 20*radius/128
%
% From Doug Noll.
% jfn added rectangular filter output
%
% $Id: fermi2d.m,v 1.4 2018/11/14 20:34:15 jfnielse Exp $

%
% [x,y]= meshdom(-64:1:63,  63:-1:-64);
% f= 1 ./ (1 + exp((sqrt(x.^2 + y.^2) - radius) ./ (10*steepness/256)));
if width < 0,
   width = 20*radius/128;
end
cent = sz/2 + 1;
x= (1:sz);
y= (1:sz)';
X= ones(size(y))*x;
Y= y*ones(size(x));

X = X-cent;
Y = Y-cent;

R = abs( X + Y.*1i );
fCirc = 1 ./ (1 + exp( (R - radius) ./ width ));

fx = 1 ./ (1 + exp( (abs(X) - radius) ./ width ));
fy = 1 ./ (1 + exp( (abs(Y) - radius) ./ width ));
fRect = fx.*fy;

return;
