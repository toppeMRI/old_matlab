function [rf, gx, gy, gz, tpad] = dispseq(nstart, nstop, varargin)
% function [rho, th, gx, gy, gz, rho1, th1, gx1, gy1, gz1, textra] = dispseq(nstart, nstop, varargin)
%
% Display pulse sequence, as specified in modules.txt, scanloop.txt, and timing.txt
%
% Inputs:
%   nstart,nstop       first and last startseq calls (as specified in scanloop.txt)
% Options:
%   'looparr'          scan loop array (see readloop.m). Default: read from scanloopfile
%   'scanloopfile'     default: 'scanloop.txt'
%   'tparams'          [start_core myrfdel daqdel timetrwait timessi]. Default: get values from timingfile.
%   'timingfile'       default: 'timing.txt'
%   'mods'             Structure containing .mod file contents (see readModules.m). Default: get values from modulesfile.
%   'modulesfile'      default: 'modules.txt'
%   'dodisplay'        true (default) or false
%
% Outputs:
%   rho                Gauss
%   th                 radians, [-pi pi]
%   gx,gy,gz           Gauss/cm

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
% (c) 2016-2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: dispseq.m,v 1.6 2018/10/09 15:22:34 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2/dispseq.m,v $

%% parse inputs

% Default values 
arg.looparr = [];
arg.scanloopfile = 'scanloop.txt';
arg.tparams = [];
arg.timingfile = 'timing.txt';
arg.mods = [];
arg.modulesfile = 'modules.txt';
arg.dodisplay = true;

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

%% read scan files as needed

% read scanloop
if isempty(arg.looparr)
	looparr = tryread(@readloop, arg.scanloopfile);
else
	looparr = arg.looparr;
end

% get timing CVs
if isempty(arg.tparams)
	TPARAMS = tryread(@readTimingFile, arg.timingfile);
else
	TPARAMS = arg.tparams;
end
[start_core myrfdel daqdel timetrwait timessi] = deal(TPARAMS(1), TPARAMS(2), TPARAMS(3), TPARAMS(4), TPARAMS(5));

% read module waveforms
if isempty(arg.mods)
	cores = tryread(@readModules, arg.modulesfile);
else
	cores = arg.mods;
end

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
	tpad = max(cores{ic}.dur - tmin, 0);                                   % silence at end of core
	tminwait = 12;   % (us) min length of wait pulse.
	if size(looparr,2)>13
		textra = looparr(it,14);
		tpad = tpad + max(textra, tminwait);    % waitcore duration (see toppev2.e)
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
	
	rho1 = [zeros(round((start_core+coredel)/dt),1); ia_rf/max_pg_iamp*  abs(cores{ic}.rf(:,waveform));  zeros(round((timetrwait+timessi)/dt),1)];
	th1  = [zeros(round((start_core+coredel)/dt),1); ia_th/max_pg_iamp*angle(cores{ic}.rf(:,waveform));  zeros(round((timetrwait+timessi)/dt),1)];
	gx1  = [zeros(round((start_core)/dt),1);         ia_gx/max_pg_iamp*gxit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gy1  = [zeros(round((start_core)/dt),1);         ia_gy/max_pg_iamp*gyit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gz1  = [zeros(round((start_core)/dt),1);         ia_gz/max_pg_iamp*gzit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];

	% apply RF phase offset
	if cores{ic}.hasRF
		th1 = th1 + looparr(it,12)/max_pg_iamp*pi;
		th1 = angle(exp(1i*th1));   % wrap to [-pi pi] range
	end

	rho = [rho; rho1; zeros(round(tpad/dt),1)];
	th  = [th;  th1;  zeros(round(tpad/dt),1)];
	gx  = [gx;  gx1;  zeros(round(tpad/dt),1)];
	gy  = [gy;  gy1;  zeros(round(tpad/dt),1)];
	gz  = [gz;  gz1;  zeros(round(tpad/dt),1)];

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
