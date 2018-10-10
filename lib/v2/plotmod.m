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

nwavs = size(b1,2);    % number of waveforms on each axis
nt = size(b1,1);
dt = 4e-3;  % ms
T = linspace(dt/2,nt*dt-dt/2,nt);

figure;
subplot(324); plot(T,rho);   ylabel('abs(rf)    G');
subplot(325); plot(T,theta); ylabel('angle(rf)  rad');
subplot(321); plot(T,gx);    ylabel('gx         G/cm');
subplot(322); plot(T,gy);    ylabel('gy         G/cm');
subplot(323); plot(T,gz);    ylabel('gz         G/cm');
xlabel('time (msec)');

return;
% EOF
