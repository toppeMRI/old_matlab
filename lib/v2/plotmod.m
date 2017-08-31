function plotmod(fname)
% function plotmod(fname)
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.

% TOPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU Library General Public License as published by
% the Free Software Foundation version 2.0 of the License.
%
% TOPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public License
% along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
% 
% (c) 2016 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod(fname);
b1 = rho.*exp(i*theta);

ncoils = size(b1,2);
nt = size(b1,1);
dt = 4e-3;  % ms
T = linspace(dt/2,nt*dt-dt/2,nt);

cols = 'bbbbbbbb';

figure;
for c = 1:ncoils
	subplot(2,ncoils,c);        plot(T,abs(b1(:,c)),cols(c));   ylabel(['coil ' num2str(c) ', abs(b1) (Gauss)']);
	xlabel('time (msec)');
	subplot(2,ncoils,c+ncoils); plot(T,angle(b1(:,c)),cols(c)); ylabel(['coil ' num2str(c) ', angle(b1)']);
	xlabel('time (msec)');
end

figure;
subplot(311); plot(T,gx); ylabel('gx (G/cm)');
subplot(312); plot(T,gy); ylabel('gy (G/cm)');
subplot(313); plot(T,gz); ylabel('gz (G/cm)');
xlabel('time (msec)');

return;
% EOF
