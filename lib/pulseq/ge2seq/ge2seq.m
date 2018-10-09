function seq = ge2seq(varargin)
% function ge2seq(modulesfile,scanloopfile,timingfile)
%
% TOPPE to Pulseq file conversion.
% Requires the Pulseq Matlab toolbox (http://pulseq.github.io/)
%
% Current verion: use matlab/lib/v2b/ functions (readmod.m, etc)
%
% Input options:
%  modulesfile       Text file listing all .mod files. Default: 'modules.txt'.
%                    The .mod files listed must exist in the Matlab working directory, i.e., 
%                    the directory from which you call this function.
%  scanloopfile      Text file specifying the MR scan loop. Default: 'scanloop.txt'
%  timingfile        Text file specifying low-level TOPPE timing parameters. Default: 'timing.txt'.
%
% Example:
%  >> ge2seq('scanloopfile', 'myloop.txt');
%
% $Id: ge2seq.m,v 1.34 2018/10/09 22:23:37 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/pulseq/pulseq-master/matlab/ge2seq/ge2seq.m,v $

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

% defaults
arg.debug = false;
arg.modulesfile  = 'modules.txt';
arg.scanloopfile = 'scanloop.txt';
arg.timingfile   = 'timing.txt';

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

debug = arg.debug;

rasterTime = 4e-6;         % TOPPE raster time (for RF, gradients, and ADC)
max_pg_iamp = 2^15-2;      % max TOPPE/GE "instruction amplitude" (signed short int)

% get TOPPE timing CVs
fid        = fopen(arg.timingfile, 'r', 'ieee-be');
s          = fgets(fid);  % skip line
s          = fscanf(fid, '%s ', 1);
start_core = fscanf(fid, '%d\n', 1);
s          = fscanf(fid, '%s ', 1);
myrfdel    = fscanf(fid, '%d\n', 1);
s          = fscanf(fid, '%s ', 1);
daqdel     = fscanf(fid, '%d\n', 1);
s          = fscanf(fid, '%s ', 1);
timetrwait = fscanf(fid, '%d\n', 1);
s          = fscanf(fid, '%s ', 1);
timessi    = fscanf(fid, '%d\n', 1);
fclose(fid);
TPARAMS = [start_core myrfdel daqdel timetrwait timessi];

% Get module information.
fid      = fopen(arg.modulesfile, 'r', 'ieee-be');
s        = fgets(fid);  % skip line
nmodules = fscanf(fid, '%d\n', 1);
s        = fgets(fid);  % skip line

for ii = 1:nmodules
	modules{ii}.fname = fscanf(fid, '%s ', 1);
	modules{ii}.dur = fscanf(fid, '%d ', 1);
	modules{ii}.hasRF = fscanf(fid, '%d ', 1);
	modules{ii}.hasDAQ = fscanf(fid, '%d\n', 1);
end
fclose(fid);

% initialize Pulseq sequence object
seq=mr.Sequence();

% load scan loop and modules
d = tryread(@readloop, arg.scanloopfile);
cores = tryread(@readModules, arg.modulesfile);

% Loop through scanloop.txt. Add each row as one Pulseq "block".
if (debug)
	nt = 20;
else
	nt = size(d,1);    % number of startseq calls
end
for ii = 1:nt
	if ~mod(ii,250)
		fprintf('.');
	end

	% get waveforms and delay for one row (one startseq call)
	[~, ~, ~, ~, rf, gxwav, gywav, gzwav, tdelay] = dispseq(ii,ii,'looparr',d,'tparams',TPARAMS,'mods',cores,'dodisplay',false);  % rf: Gauss; gradients: Gauss/cm; tdelay: microsec

	gx = makeArbitraryGrad('x',g2pulseq(gxwav,rasterTime,seq));
	gy = makeArbitraryGrad('y',g2pulseq(gywav,rasterTime,seq));
	gz = makeArbitraryGrad('z',g2pulseq(gzwav,rasterTime,seq));

	freqOffset  = d(ii,15);                         % Hz

	module = modules{d(ii,1)};
	
	if module.hasRF
		phaseOffset = d(ii,12)/max_pg_iamp*pi;          % radians
		flip = d(ii,2)/max_pg_iamp*pi;
		rf = makeArbitraryRf(rf2pulseq(rf,rasterTime,seq), flip, 'FreqOffset', freqOffset, 'PhaseOffset', phaseOffset);
		seq.addBlock(rf,gx,gy,gz);
		if debug
			subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal))));
			subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal))));
		end
	elseif module.hasDAQ
		phaseOffset = d(ii,13)/max_pg_iamp*pi;          % radians
		adc = makeAdc(numel(gxwav),'Dwell',rasterTime,'delay',0,...
				'freqOffset',freqOffset,'phaseOffset',phaseOffset);
		seq.addBlock(gx,gy,gz,adc);
	else
		seq.addBlock(gx,gy,gz);
	end

	if debug
		subplot(2,2,[3 4]); plot(gx.waveform,'r'); hold on; plot(gy.waveform,'g'); plot(gz.waveform,'b'); hold off; title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		%input('press any key to continue');
	end

	% add delay block
	if tdelay > 0
		del = makeDelay(tdelay*1e-6); 
		seq.addBlock(del);
	end
	
end
fprintf('\n');

seq.plot();
seq.write('out.seq');

return;
	
function g = g2pulseq(g,rasterTime,seq)
% convert gradient from Gauss/cm to Hz/m, and interpolate to GradRasterTime
gamma = 42.576e6;      % Hz/T
g = g * gamma / 1.0e2;   % Hz/m
n = numel(g);
gradRasterTime = seq.gradRasterTime;
g = interp1(1:n,g,floor(gradRasterTime/rasterTime):gradRasterTime/rasterTime:n);
return;

function rf = rf2pulseq(rf,rasterTime,seq)
% convert rf units from Gauss to Hz, and interpolate to rfRasterTime
gamma = 4.2576e3;       % Hz/G
rf = rf*gamma;          % Hz
%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/rfRasterTime,L,cutoff);      % upsample from 4us to 1us
n = numel(rf);
rfRasterTime = seq.rfRasterTime;
rf = interp1(1:n,rf,floor(rfRasterTime/rasterTime):rfRasterTime/rasterTime:(n-rfRasterTime/rasterTime),'linear','extrap');
return;

