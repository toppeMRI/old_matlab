function ge2seq(modulesfile,scanloopfile,timingfile,debug)
% function ge2seq(modulesfile,scanloopfile,timingfile,[debug])
%
% TOPPE to Pulseq file conversion.
% Requires the Pulseq Matlab toolbox (http://pulseq.github.io/)
%
% Inputs:
%  modulesfile       Text file listing all .mod files (modules.txt).
%                    The .mod files listed must exist in the Matlab working directory, i.e., 
%                    the directory from which you call this function.
%  scanloopfile      Text file specifying the MR experiment (scanloop.txt)
%  timingfile        Text file specifying low-level TOPPE timing parameters (timing.txt).
%
% Example:
%  >> ge2seq('modules.txt','scanloop.txt','timing.txt')

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
% $Id: ge2seq.m,v 1.16 2017/09/15 19:53:29 jfnielse Exp $

if ~exist('debug','var')
	debug = false;
end

rasterTime = 4e-6;         % TOPPE raster time (for RF, gradients, and ADC)
max_pg_iamp = 2^15-2;      % max TOPPE/GE "instruction amplitude" (signed short int)

% get TOPPE timing CVs
fid = fopen(timingfile, 'r', 'ieee-be');
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
TPARAMS = [start_core myrfdel daqdel timetrwait timessi];

% Get module information.
fid = fopen(modulesfile, 'r', 'ieee-be');
s = fgets(fid);  % skip line
nmodules = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
ishp = 1;  % index into shapeLibrary
for ii = 1:nmodules
	modules{ii}.fname = fscanf(fid, '%s ', 1);
	modules{ii}.dur = fscanf(fid, '%d ', 1);
	modules{ii}.hasRF = fscanf(fid, '%d ', 1);
	modules{ii}.hasDAQ = fscanf(fid, '%d\n', 1);
	%[desc rho th gx gy gz]	= readmod(modules{ii}.fname,false);
end
fclose(fid);

% initialize Pulseq sequence object
seq=mr.Sequence();

% Loop through scanloop.txt. Add each row as one Pulseq "block".
d = readloop(scanloopfile);
nt = size(d,1);    % number of startseq calls
%for ii = 1:nt   %82
for ii = 1:100
	if ~mod(ii,2*round(nt/10/2))
		fprintf('.');
	end
	% get waveforms and delay for one row (one startseq call)
	[~,~,~,~,~,rho,th,gxwav,gywav,gzwav,delay] = scansim(ii,ii,d,TPARAMS,false);  % rf: Gauss; gradients: Gauss/cm; delay: microsec

	gx = makeArbitraryGrad('x',g2pulseq(gxwav,rasterTime));
	gy = makeArbitraryGrad('y',g2pulseq(gywav,rasterTime));
	gz = makeArbitraryGrad('z',g2pulseq(gzwav,rasterTime));

	freqOffset  = d(ii,15);                         % Hz
	phaseOffset = d(ii,12)/max_pg_iamp*pi;          % radians

	module = modules{d(ii,1)};
	
	if module.hasRF
		flip = d(ii,2)/max_pg_iamp*pi;
		rf = makeArbitraryRf(rf2pulseq(rho.*exp(1i*th),rasterTime),flip);
		seq.addBlock(rf,gx,gy,gz);
		if debug
			subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal))));
			subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal))));
			%subplot(221); plot(rho,'r');
			%subplot(222); plot(th);
		end
	elseif module.hasDAQ
		adc = makeAdc(numel(gxwav),'Dwell',rasterTime,'delay',0,...
							'freqOffset',freqOffset, 'phaseOffset',phaseOffset);
		seq.addBlock(gx,gy,gz,adc);
	else
		seq.addBlock(gx,gy,gz);
	end

	if debug
		subplot(2,2,[3 4]); plot(gx.waveform,'r'); hold on; plot(gy.waveform,'g'); plot(gz.waveform,'b'); hold off; title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		input('press any key to continue');
	end

	% add delay block
	if delay > 0
		del = makeDelay(delay*1e-6); 
		seq.addBlock(del);
	end
	
end
fprintf('\n');

seq.plot();
seq.write('test.seq');

return;
	
function g = g2pulseq(g,rasterTime)
% convert gradient Gauss/cm to Hz/cm, and interpolate to GradRasterTime
gamma = 42.576e6;       % Hz/T
g = 1e-2 * g * gamma;   % Gauss/cm
GradRasterTime = 1e-5;  % Pulseq gradient raster time
n = numel(g);
g = interp1(1:n,g,floor(GradRasterTime/rasterTime):GradRasterTime/rasterTime:n);
return;

function rf = rf2pulseq(rf,rasterTime)
% convert rf units from Gauss to Hz, and interpolate to RfRasterTime
gamma = 42.576e6;   % Hz/T
rf = 1e-4*rf*gamma;
RfRasterTime = 1e-6;    % Pulseq RF raster time
%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/RfRasterTime,L,cutoff);      % upsample from 4us to 1us
n = numel(rf);
rf = interp1(1:n,rf,floor(RfRasterTime/rasterTime):RfRasterTime/rasterTime:(n-RfRasterTime/rasterTime),'linear','extrap');
return;

