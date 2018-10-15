function seq = ge2seq(varargin)
% function ge2seq(modulesfile,scanloopfile,timingfile)
%
% TOPPE to Pulseq file conversion.
% Requires the Pulseq Matlab toolbox (http://pulseq.github.io/)
%
% Current TOPPE functions: matlab/lib/v2b/ functions (readmod.m, etc)
%
% Input options:
%  modulesfile     Text file listing all .mod files. Default: 'modules.txt'.
%                  The .mod files listed must exist in the Matlab path.
%  scanloopfile    Text file specifying the MR scan loop. Default: 'scanloop.txt'
%  timingfile      Text file specifying low-level TOPPE timing parameters. Default: 'timing.txt'.
%
% Examples:
%  >> ge2seq();
%  >> ge2seq('debug', 'true');
%  >> ge2seq('scanloopfile', 'myloop.txt');
%
% $Id: ge2seq.m,v 1.47 2018/10/15 11:07:02 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/pulseq/pulseq-master/matlab/ge2seq/ge2seq.m,v $

% defaults
arg.debug = false;
arg.debugAdc = false;
arg.modulesfile  = 'modules.txt';
arg.scanloopfile = 'scanloop.txt';
arg.timingfile   = 'timing.txt';

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

%% initialize Pulseq sequence object
lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
               'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
               'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  
seq = mr.Sequence(lims);

%% TOPPE scan info
fid      = fopen(arg.modulesfile, 'r', 'ieee-be');
s        = fgets(fid);  % skip line
nmodules = fscanf(fid, '%d\n', 1);
s        = fgets(fid);  % skip line
for ii = 1:nmodules
	modules{ii}.fname  = fscanf(fid, '%s ', 1);
	modules{ii}.dur    = fscanf(fid, '%d ', 1);
	modules{ii}.hasRF  = fscanf(fid, '%d ', 1);
	modules{ii}.hasDAQ = fscanf(fid, '%d\n', 1);
end
fclose(fid);

geRasterTime = 4e-6;        % TOPPE raster time (for RF, gradients, and ADC)
max_pg_iamp  = 2^15-2;      % max TOPPE/GE "instruction amplitude" (signed short int)

d       = tryread(@readloop,       arg.scanloopfile);    % scanloop array
cores   = tryread(@readModules,    arg.modulesfile);     % module waveforms
TPARAMS = tryread(@readTimingFile, arg.timingfile);      % low-level (fixed) timing parameters

%% Loop through scanloop.txt. Add each row as one Pulseq "block".
nt = size(d,1);    % number of startseq calls
for ii = 1:nt
	if ~mod(ii,250)
		fprintf('.');
	end

	% get waveforms and delay for one row (one startseq call)
	% rf: Gauss; gradients: Gauss/cm; tdelay: microsec
	[~, ~, ~, ~, rfwav, gxwav, gywav, gzwav, tdelay] = dispseq(ii,ii,'looparr',d,...
		'tparams',TPARAMS,'mods',cores,'dodisplay',false);  

	% pulseq likes row vectors
	rfwav = rfwav(:)';
	gxwav = gxwav(:)';
	gywav = gywav(:)';
	gzwav = gzwav(:)';

	% convert to Pulseq units and rastertimes
	% rf:   Hz,   1us
   % grad: Hz/m, 10us
	rfwavPulseq = rf2pulseq(rfwav,geRasterTime,seq);
	gxwavPulseq = g2pulseq( gxwav,geRasterTime,seq);
	gywavPulseq = g2pulseq( gywav,geRasterTime,seq);
	gzwavPulseq = g2pulseq( gzwav,geRasterTime,seq);

	% ensure equal duration (interpolation to Pulseq rastertimes can result in unequal duration)
	% not needed?
	trf   = length(rfwavPulseq) * seq.rfRasterTime;
	tgrad = length(gxwavPulseq) * seq.gradRasterTime;
	ngradextra = ceil((trf-tgrad)/seq.gradRasterTime);
	gxwavPulseq = makeevenlength( [gxwavPulseq zeros(1, ngradextra)] );
	gywavPulseq = makeevenlength( [gywavPulseq zeros(1, ngradextra)] );
	gzwavPulseq = makeevenlength( [gzwavPulseq zeros(1, ngradextra)] );
	tgrad  = length(gxwavPulseq) * seq.gradRasterTime;
	rfwavPulseq = [rfwavPulseq zeros(1,round((tgrad-trf)/seq.rfRasterTime))];

	% Make Pulseq gradient structs (even all zero waveforms)
	gx = mr.makeArbitraryGrad('x', gxwavPulseq, lims);
	gy = mr.makeArbitraryGrad('y', gywavPulseq, lims);
	gz = mr.makeArbitraryGrad('z', gzwavPulseq, lims);

	% bitmask indicating non-zero gradients
	hasg = 0;   
	if ~all(gxwavPulseq == 0)
		hasg = bitset(hasg,1);
	end
	if ~all(gywavPulseq == 0)
		hasg = bitset(hasg,2);
	end
	if ~all(gzwavPulseq == 0)
		hasg = bitset(hasg,3);
	end
	strArg = getArgStr(hasg);        % 'gz' or 'gx,gy,gz' or... as appropriate

	freqOffset  = d(ii,15);                         % Hz
	module = modules{d(ii,1)};

	if module.hasRF
		phaseOffset = d(ii,12)/max_pg_iamp*pi;          % radians
		flip = pi; % d(ii,2)/max_pg_iamp*pi;
		rf = mr.makeArbitraryRf(rfwavPulseq, flip, 'FreqOffset', freqOffset, ...
			'PhaseOffset', phaseOffset, 'system', lims);

		if isempty(strArg)
			seq.addBlock(rf);
		else
			eval( sprintf( 'seq.addBlock(rf, %s)', strArg) );
		end

		if arg.debug
			clf;
			subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal)))); ylabel('Hz');
			subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal)))); ylabel('rad');
		end
	elseif module.hasDAQ
		phaseOffset = d(ii,13)/max_pg_iamp*pi;          % radians
		if arg.debugAdc
			% delay and shorten adc window
			nadc = 2*round(numel(gx.waveform)/4);
			delay = round(nadc/4)*seq.gradRasterTime;
			adc = mr.makeAdc(nadc, lims, 'Dwell', seq.gradRasterTime, 'delay', delay, ...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		else
			adc = mr.makeAdc(numel(gx.waveform), lims, 'Dwell', seq.gradRasterTime, 'delay',0,...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		end

		if isempty(strArg)
			seq.addBlock(adc);
		else
			eval( sprintf( 'seq.addBlock(%s, adc)', strArg) );
		end
	else
		if ~isempty(strArg)
			eval( sprintf( 'seq.addBlock(%s)', strArg) );
		end
	end

	if arg.debug
		if ~module.hasRF
			clf;
		end
		subplot(2,2,[3 4]); 
		hold on
		if bitget(hasg, 1)
			plot(gx.waveform,'r'); ylabel('Hz/m'); hold on; 
		end
		if bitget(hasg, 2)
			plot(gy.waveform,'g'); hold on;
		end
		if bitget(hasg, 3)
			plot(gz.waveform,'b'); 
		end
		if ~hasg
			clf(sfh);
		end
			
		hold off; %title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		input('press any key to continue');
	end

	% add delay block
	if tdelay > 0
		del = mr.makeDelay(round(tdelay*1e-6,5)); %delay also needs to be in multiples of rastertimes of 10us
		seq.addBlock(del);
	end
	
end
fprintf('\n');

%seq.plot();
seq.write('out.seq');

return;

	
%% convert gradient from Gauss/cm to Hz/m, and interpolate to seq.gradRasterTime
function gout = g2pulseq(g,geRasterTime,seq)
gamma = 4.2576e3;      % Hz/G
g = g * gamma * 100;   % Hz/m
T = numel(g)*geRasterTime;    % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.gradRasterTime:(T-seq.gradRasterTime);
gout = interp1(tge, g, t, 'linear', 'extrap');
return;

%% convert rf units from Gauss to Hz, and interpolate to seq.rfRasterTime
function rfout = rf2pulseq(rf,geRasterTime,seq)
gamma = 4.2576e3;       % Hz/G
rf = rf*gamma;          % Hz
T = numel(rf)*geRasterTime;   % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.rfRasterTime:(T-seq.rfRasterTime);
rfout = interp1(tge, rf, t, 'linear', 'extrap');
%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/rfRasterTime,L,cutoff);      % upsample from 4us to 1us
return;

%% get gradient arguments (as string) to pass to seq.addBlock()
function argStr = getArgStr(hasg)

switch hasg
	case 0
		argStr = '';  % no gradients
	case 1
		argStr = 'gx'; 
	case 2
		argStr = 'gy'; 
	case 4
		argStr = 'gz'; 
	case 3
		argStr = 'gx, gy'; 
	case 5
		argStr = 'gx, gz'; 
	case 6
		argStr = 'gy, gz'; 
	case 7
		argStr = 'gx, gy, gz'; 
end

return;
