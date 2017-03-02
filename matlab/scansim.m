function [rho, th, gx, gy, gz] = scansim(nstart, nstop, d)
% function [rho, th, gx, gy, gz] = scansim(nstart, nstop, [d])
% display pulse sequence, as specified in cores.txt, scanloop.txt, and timing.txt
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
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
%
% $Id: scansim.m,v 1.27 2017/02/07 19:44:15 jfnielse Exp $

if ~exist('nstart','var')
	nstart = 1;
end
if ~exist('nstop','var')
	nstop = nstart+10;
end

% get timing CVs
fid = fopen('timing.txt', 'r', 'ieee-be');
s = fgets(fid);  % skip line
s = fscanf(fid, '%s ', 1);
start_core = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
myrfdel = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
daqdel = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
timetrwait = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
timessi = fscanf(fid, '%d\n', 1);
fclose(fid);

% get waveforms
fid = fopen('cores.txt', 'r', 'ieee-be');
s = fgets(fid);  % skip line
ncores = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
for ic = 1:ncores
	cores{ic}.fname = fscanf(fid, '%s ', 1);
	cores{ic}.dur = fscanf(fid, '%d ', 1);
	cores{ic}.hasRF = fscanf(fid, '%d ', 1);
	cores{ic}.hasDAQ = fscanf(fid, '%d\n', 1);
	[desc,cores{ic}.rho,cores{ic}.th,cores{ic}.gx,cores{ic}.gy,cores{ic}.gz,cores{ic}.paramsint16,cores{ic}.paramsfloat] ...
		= readmod(cores{ic}.fname,false);
	cores{ic}.wavdur = numel(cores{ic}.gx)*4;   % waveform duration [us]
end
fclose(fid);

% get scan loop
if ~exist('d', 'var')
	d = readloop('scanloop.txt');
end

% build sequence. each sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
dt = 4;  % us
max_pg_iamp = 2^15-2;
for it = nstart:nstop
	ic = d(it,1);   % core id 
	if cores{ic}.hasRF
		ia_rf = d(it,2);
	else
		ia_rf = 0;
	end
	ia_th = d(it,3);
	ia_gx = d(it,4);
	ia_gy = d(it,5);
	ia_gz = d(it,6);

	if cores{ic}.hasRF
		coredel = myrfdel;
	elseif cores{ic}.hasDAQ
		coredel = daqdel;
	else
		coredel = 0;
	end

	tmin = start_core + coredel + cores{ic}.wavdur + timetrwait + 8 + timessi;   % mimimum core duration. 8 is for wait pulse.
	textra = max(cores{ic}.dur - tmin, 0);                                       % silence at end of core
	if size(d,2)>13
		textra = textra + d(it,14);
	end

	% apply in-plane (xy) rotation
	iphi = d(it,11);
	phi = iphi/max_pg_iamp*pi;    % rad, [-pi pi]
	Gxy = [cos(phi) -sin(phi); sin(phi) cos(phi)]*[cores{ic}.gx(:)'; cores{ic}.gy(:)'];
	gxit = Gxy(1,:)';
	gyit = Gxy(2,:)';

	rho = [rho; zeros(round((start_core+coredel)/dt),1); ia_rf/max_pg_iamp*cores{ic}.rho; zeros(round((timetrwait+timessi)/dt),1);         zeros(round(textra/dt),1)];
	th  = [th;  zeros(round((start_core+coredel)/dt),1); ia_th/max_pg_iamp*cores{ic}.th;  zeros(round((timetrwait+timessi)/dt),1);         zeros(round(textra/dt),1)];
	gx  = [gx;  zeros(round((start_core)/dt),1);         ia_gx/max_pg_iamp*gxit;          zeros(round((timetrwait+timessi+coredel)/dt),1); zeros(round(textra/dt),1)];
	gy  = [gy;  zeros(round((start_core)/dt),1);         ia_gy/max_pg_iamp*gyit;          zeros(round((timetrwait+timessi+coredel)/dt),1); zeros(round(textra/dt),1)];
	gz  = [gz;  zeros(round((start_core)/dt),1);         ia_gz/max_pg_iamp*cores{ic}.gz;  zeros(round((timetrwait+timessi+coredel)/dt),1); zeros(round(textra/dt),1)];

	%fprintf(1, 'it %d: tmin = %.3f ms, rf t = %.3f ms, grad t = %.3f ms\n', it, tmin/1000, numel(rho)*dt*1e-3, numel(gx)*dt*1e-3);
end

% plot
T = (0:(numel(rho)-1))*dt/1000; % msec
gmax = 5;  % Gauss/cm
subplot(511); plot(T, rho); ylabel('rho');   axis([T(1) 1.01*T(end) -1.1*min(rho) 1.1*max(rho)]);
subplot(512); plot(T, th);  ylabel('theta'); axis([T(1) 1.01*T(end) -1.3*pi 1.3*pi]);
subplot(513); plot(T, gx);  ylabel('gx'); axis([T(1) 1.01*T(end) -1.2*gmax 1.2*gmax]);;
subplot(514); plot(T, gy);  ylabel('gy'); axis([T(1) 1.01*T(end) -1.2*gmax 1.2*gmax]);;
subplot(515); plot(T, gz);  ylabel('gz'); axis([T(1) 1.01*T(end) -1.2*gmax 1.2*gmax]);;
xlabel('msec');

return;


% EOF
