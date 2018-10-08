function [looparr, mods, tparams] = readAllFiles(varargin)
% function [looparr, mods, tparams] = readAllFiles(varargin)
%
% Display pulse sequence, as specified in modules.txt, scanloop.txt, and timing.txt
%
% Input options:
% Options:
%   'scanloopfile'     default: 'scanloop.txt'
%   'timingfile'       default: 'timing.txt'
%   'modulesfile'      default: 'modules.txt'
%
% Outputs:
%   looparr            [nt 16] array containing scan loop parameters (see readloop.m)
%   mods               struct containing .mod file contents (see readModules.m)
%   tparams            [1 5] array containing low-level sequence timing parameters (see readTimingFile.m)

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
% (c) 2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: readAllFiles.m,v 1.1 2018/10/08 14:14:01 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2/readAllFiles.m,v $

%% parse inputs

% Default values 
arg.scanloopfile = 'scanloop.txt';
arg.timingfile = 'timing.txt';
arg.modulesfile = 'modules.txt';

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

looparr = tryread(@readloop, arg.scanloopfile);

% low-level timing CVs
TPARAMS = tryread(@readTimingFile, arg.timingfile);

% read module waveforms
mods = tryread(@readModules, arg.modulesfile);

% build sequence. each sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
dt = 4;  % us
max_pg_iamp = 2^15-2;
for it = nstart:nstop
	ic = looparr(it,1);   % core id 
	if cores{ic}.hasRF
		ia_rf = looparr(it,2);
	else
		ia_rf = 0;
	end
	ia_th = looparr(it,3);
	ia_gx = looparr(it,4);
	ia_gy = looparr(it,5);
	ia_gz = looparr(it,6);

	if cores{ic}.hasRF
		coredel = myrfdel;
	elseif cores{ic}.hasDAQ
		coredel = daqdel;
	else
		coredel = 0;
	end

	tmin = start_core + coredel + cores{ic}.wavdur + timetrwait + timessi;   % mimimum core duration (us). 
	textra = max(cores{ic}.dur - tmin, 0);                                   % silence at end of core
	tminwait = 12;   % (us) min length of wait pulse.
	if size(looparr,2)>13
		textra = textra + max(looparr(it,14),tminwait);    % waitcore duration (see toppev2.e)
	end

	waveform = looparr(it,16);

	% get gradients and apply in-plane (xy) rotation
	gxit = cores{ic}.gx(:,waveform);
	gyit = cores{ic}.gy(:,waveform);
	gzit = cores{ic}.gz(:,waveform);
	iphi = looparr(it,11);
	phi = iphi/max_pg_iamp*pi;    % rad, [-pi pi]
	Gxy = [cos(phi) -sin(phi); sin(phi) cos(phi)]*[gxit(:)'; gyit(:)'];
	gxit = Gxy(1,:)';
	gyit = Gxy(2,:)';
	
	rho1 = [zeros(round((start_core+coredel)/dt),1); ia_rf/max_pg_iamp*cores{ic}.rho(:,waveform); zeros(round((timetrwait+timessi)/dt),1)];
	th1  = [zeros(round((start_core+coredel)/dt),1); ia_th/max_pg_iamp*cores{ic}.th(:,waveform);  zeros(round((timetrwait+timessi)/dt),1)];
	gx1  = [zeros(round((start_core)/dt),1);         ia_gx/max_pg_iamp*gxit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gy1  = [zeros(round((start_core)/dt),1);         ia_gy/max_pg_iamp*gyit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gz1  = [zeros(round((start_core)/dt),1);         ia_gz/max_pg_iamp*gzit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];

	% apply RF phase offset
	if cores{ic}.hasRF
		th1 = th1 + looparr(it,12)/max_pg_iamp*pi;
		th1 = angle(exp(1i*th1));   % wrap to [-pi pi] range
	end

	rho = [rho; rho1; zeros(round(textra/dt),1)];
	th  = [th;  th1;  zeros(round(textra/dt),1)];
	gx  = [gx;  gx1;  zeros(round(textra/dt),1)];
	gy  = [gy;  gy1;  zeros(round(textra/dt),1)];
	gz  = [gz;  gz1;  zeros(round(textra/dt),1)];

	%fprintf(1, 'it %d: tmin = %.3f ms, rf t = %.3f ms, grad t = %.3f ms\n', it, tmin/1000, numel(rho)*dt*1e-3, numel(gx)*dt*1e-3);
end

% plot
if arg.dodisplay
	T = (0:(numel(rho)-1))*dt/1000; % msec
	gmax = 5;  % Gauss/cm
    srho = max(1.1*max(abs(rho(:))),0.05);
	subplot(511); plot(T, rho); ylabel('rho');   axis([T(1) 1.01*T(end) -srho srho]);
	subplot(512); plot(T, th);  ylabel('theta'); axis([T(1) 1.01*T(end) -1.3*pi 1.3*pi]);
	subplot(513); plot(T, gx);  ylabel('gx'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	%gmax = 1;  % Gauss/cm
	subplot(514); plot(T, gy);  ylabel('gy'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	subplot(515); plot(T, gz);  ylabel('gz'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	xlabel('msec');
end

return;


% EOF
