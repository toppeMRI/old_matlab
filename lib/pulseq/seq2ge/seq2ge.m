function seq2ge(seqfile)
% function seq2ge(seqfile)
%
% Convert Pulseq .seq file to TOPPE files, for playing sequence on GE scanner.
%
% The TOPPE pulse sequence needs the following files:
%   *.mod:          one .mod file corresponds to one "unique" block (see below)
%   modules.txt:    list of .mod files, and flags indicating whether the .mod file is RF/ADC/(gradients only)
%   scanloop.txt:   sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)

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
% $Id: loop2txt.m,v 1.1 2017/08/29 16:00:17 jfnielse Exp $

dotroubleshoot = false;

% read .seq file
seq=mr.Sequence();
seq.read(seqfile);

%mxslew = 15;  % gradient slew rate limit.

if 1

%
%
% Find blocks that are unique in terms of waveforms and timing (i.e., waveform amplitudes, RF/ADC phase, etc can differ).
% The resulting set of unique blocks will be referred to as "cores", to distinguish them from the blocks defined in the .seq file.
% For each core, one .mod file is generated.
%
%

fprintf(1, '%s: Finding unique cores... ', seqfile);

cores = [1 seq.blockEvents(1,:)];   % [block# size(seq.blockEvents,2)] Array containing unique blocks (cores)
%coreinds = 1:size(seq.blockEvents,1);   

% loop through blocks
for ib= 1:size(seq.blockEvents,1)

	% Check to see if block ib is unique with respect to previously defined cores.

	Isunique = 0; % counts occurrences of uniqueness w.r.t. previously identified cores

	for ic = 1:size(cores,1)
		isunique = false;
		Idiff = find(seq.blockEvents(ib,:) ~= cores(ic,2:end)) ;
		if ~isempty(Idiff)
			% Found block with at least one different event ID.
			% Now we check to see if at least one of the events is unique.
			block = seq.getBlock(ib); 
			core = seq.getBlock(cores(ic,1));   % block already in cores array

			% TODO: check uniqueness of RF events
			% Find unique RF events
			%nrfevents = 0;
			%for ii=1:seq.rfLibrary.Count
			%end
			% For now, assume each RF event is unique: (TODO: fix this)
			if seq.blockEvents(ib,2) ~= cores(ic,3)
				isunique = true;
			end

			% each delay event is unique (by definition)
			if seq.blockEvents(ib,1) ~= cores(ic,2)
				isunique = true;
			end

			% check uniqueness of gradient events
			gradChannels={'gx','gy','gz'};
			for j=1:length(gradChannels)
				gradb = block.(gradChannels{j});
				gradc = core.(gradChannels{j});
				if xor(~isempty(gradb),~isempty(gradc))
					isunique = true;  % count as unique
				elseif ~isempty(gradb)
					if ~strcmp(gradb.type,gradc.type)
						% gradients are different type, so assumed to be unique
						isunique = true;
					else
						if strcmp(gradb.type,'trap')
							if any([gradb.riseTime gradb.flatTime gradb.fallTime] ~= [gradc.riseTime gradc.flatTime gradc.fallTime])
								isunique = true; % count as unique
							end
						else
							% TODO: deal with non-trap gradients here
						end
					end
				end
			end
		end

		Isunique(ic) = isunique; 
	end

	% coreinds = mapping from block # to corresponding (unique) core #
	if all(Isunique) 
		% Block ib is unique w.r.t. all previously known cores, so add block ib to list of cores.
		cores = [cores; [ib seq.blockEvents(ib,:)]];
		coreinds(ib) = size(cores,1);
	else
		% Block is not unique. 
		isim = find(~Isunique);
		coreinds(ib) = isim;
	end

end

fprintf(1, 'done.\n');
		
save coreinfo cores coreinds

else

load coreinfo   % for development

end

fprintf(1,'Found the following unique cores: \n\n   block# delay  RF    gx    gy    gz    ADC\n');
disp(cores);

% find maximum gradient amplitudes for each core
fprintf(1, 'Finding maximum gradient amplitudes for each core... ');
for ic = 1:size(cores,1)
	gxmax = 0;
	gymax = 0;
	gzmax = 0;

	% loop through all blocks that "belong" to core ic
	inds = find(coreinds == ic);
	for ii = inds  
		b = seq.getBlock(ii);
		if ~isempty(b.gx)
			if abs(b.gx.amplitude) > gxmax
				gxmax = abs(b.gx.amplitude);   % Hz/cm
			end
		end
		if ~isempty(b.gy)
			if abs(b.gy.amplitude) > gymax
				gymax = abs(b.gy.amplitude);
			end
		end
		if ~isempty(b.gz)
			if abs(b.gz.amplitude) > gzmax
				gzmax = abs(b.gz.amplitude);
			end
		end
	end

	Gxmax(ic) = g2gcm(gxmax);   % Gauss/cm
	Gymax(ic) = g2gcm(gymax);
	Gzmax(ic) = g2gcm(gzmax);
end

% For troubleshooting, reduce gradient amplitudes.
if dotroubleshoot
redfactor = 1/2;
Gxmax(3:4) = Gxmax(3:4)*redfactor;
Gymax(3:4) = Gymax(3:4)*redfactor;
Gzmax(3:4) = Gzmax(3:4)*redfactor;
end

fprintf(1, 'done.\n');
fprintf(1,'Found the following peak gradient (Gauss/cm) for each core: \n\n       core1     core2    ...\n');
disp([Gxmax ; Gymax; Gzmax]);





%
%
% Write each core to a .mod file, and create modules.txt.
%
%

fprintf(1, 'Writing .mod files and modules.txt... ');

fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n',size(cores,1));
fprintf(fid,'wavfile_name    duration (us)     has_RF?     has_ADC?\n');

% loop through cores
for ic = 1:size(cores,1)

	core = seq.getBlock(cores(ic,1)); 

	% defaults
	hasrf = 0;
	hasadc= 0;
	rf = [];
	gx = [];
	gy = [];
	gz = [];

	dt = 4e-6;   % dwell time (GE)

	% rf waveform
	if ~isempty(core.rf)
		hasrf = 1;
		rf = decimate(core.rf.signal,round(dt/1e-6));  % downsample from 1us to dwell time
		gamma = 4.2576e7;   % Hz/T
		rf = 1e4*rf/gamma; % Gauss
	end

	% data acquisition?
	if ~isempty(core.adc)
		hasadc = 1;
	end

	if hasrf & hasadc
		error('Cannot transmit RF and acquire data in same block. Redesign the .seq file.');
	end
	
	% gradients
	gradChannels={'gx','gy','gz'};
	gwav = [];
	for j=1:length(gradChannels)
		gc  = core.(gradChannels{j});
		switch gradChannels{j}
			case 'gx'
				gmax = Gxmax(ic);  % Gauss/cm
			case 'gy'
				gmax = Gymax(ic);
			case 'gz'
				gmax = Gzmax(ic);
		end
		if ~isempty(gc)
			if strcmp(gc.type,'trap')
				gwav = [linspace(0,gmax,ceil(gc.riseTime/dt))  gmax*ones(1,floor(gc.flatTime/dt))  linspace(gmax,0,ceil(gc.fallTime/dt))];
			else
				% TODO: deal with non-trap gradient here
			end
		end
		switch gradChannels{j}
			case 'gx'
				gx = gwav;  % Gauss/cm
			case 'gy'
				gy = gwav;
			case 'gz'
				gz = gwav;
		end
	end

	% number of waveform samples must be the same for all waveforms
	nsamples = max([numel(rf) numel(gx) numel(gy) numel(gz)]);
	if nsamples == 0  % delay core
		% put dummy pulse in .wav file (needed for RF parameter calculations)
		ndummy=100;   
		rf = 0.01*ones(ndummy,1);  % RF waveform can't be all zeroes
		gx = 0*ones(ndummy,1);
		gy = gx;
		gz = gy;
	else
		rf = [rf(:); zeros(nsamples-numel(rf),1)];
		gx = [gx(:); zeros(nsamples-numel(gx),1)];
		gy = [gy(:); zeros(nsamples-numel(gy),1)];
		gz = [gz(:); zeros(nsamples-numel(gz),1)];
	end

	% RF waveform can't be all zeroes (this has to do with the fact that both rf and gradients-only blocks use the .mod file format)
	if isempty(core.rf)
		rf = 0.01*ones(size(rf));  
	end

	% Done creating rf, gx, gy, gz.
	% Now write core to .wav file.

	if ~isempty(core.adc)
		fname = 'readout.mod';
	elseif ~isempty(core.rf)
		fname = 'tipdown.mod';
	else
		fname = sprintf('core%d.mod', ic);
	end
	nomflip = 90;
	if nsamples == 0
		type = 'delay';
	elseif ~isempty(core.rf)
		type = 'excitation';
	elseif ~isempty(core.adc)
		type = 'acquisition';
	else
		type = 'gradients only';
	end

	rf = makeGElength(rf);  % make sure waveform length is divisible by 4 (pad with zeroes at end if necessary)
	gx = makeGElength(gx); 
	gy = makeGElength(gy);
	gz = makeGElength(gz);

	desc = sprintf('Created with seq2ge. Type: %s.', type);
	hdrfloats = [];
	if ~isempty(core.adc)
		if rem(core.adc.dwell,dt)
			error('dwell time must be multiple of GE sampling rate (4us)');
		else
			decimation = core.adc.dwell/dt;  % better be integer
		end
		npre = round(core.adc.delay/dt);  % number of 4us gradient samples before start of ADC
		nadc = core.adc.numSamples*decimation;   % number of samples that the GE sequence acquires (dwell > 4us is handled in recon)
		hdrints = [0 numel(rf) npre nadc]
		hdrints(10) = decimation;
	else
		hdrints = [0 numel(rf)];
	end
	addrframp = false;
	
	%mat2wav(abs(rf),angle(rf),gx,gy,gz,nomflip,fname,desc,hdrfloats,hdrints,addrframp);
	mat2mod(abs(rf),angle(rf),gx,gy,gz,nomflip,fname,desc,hdrfloats,hdrints,addrframp);

	% write entry in modules.txt for this core
	if isempty(core.delay)
		duration = round(dt*numel(rf)*1e6); % us
	else
		duration = dt*1e6*round(core.delay.delay/dt);  % us (must be multiple of dwell time)
	end
	fprintf(fid,'%s\t%d\t%d\t%d\n', fname, duration, hasrf, hasadc);	
end

fclose(fid);

fprintf(1, 'done.\n');



%
%
% Write scanloop.txt, which specifices the scan sequence (along with modules.txt and the .mod files).
%
%

fprintf(1, 'Writing scanloop.txt... ');

dabon = 1;
daboff = 0;
max_pg_iamp = 2^15-2;  
ia_th = max_pg_iamp;
phi = 0;  % in-plane gradient rotation. Useful for, e.g., multishot spiral.

dabecho = 0;

ia_rf = max_pg_iamp;

dabview = 0;

d = []; 

% loop through blocks
for ib = 1:size(seq.blockEvents,1)

	b = seq.getBlock(ib); 

	ic = coreinds(ib);

	irfphase = 0;   % for now. TODO
	
	% set waveform amplitudes
	if ~isempty(b.gx)
		ia_gx = 2*round(max_pg_iamp/2 * g2gcm(b.gx.amplitude) / Gxmax(ic) );
	else
		ia_gx = 0;
	end

	if ~isempty(b.gy)
		ia_gy = 2*round(max_pg_iamp * g2gcm(b.gy.amplitude) / Gymax(ic) /2 );
	else
		ia_gy = 0;
	end

	if ~isempty(b.gz)
		ia_gz = 2*round(max_pg_iamp * g2gcm(b.gz.amplitude) / Gzmax(ic) /2 );
	else
		ia_gz = 0;
	end

	if ~isempty(b.adc)
		dabmode = dabon;
		dabview = dabview + 1; 	% For now, store data in increasing view order. TODO
	else
		dabmode = daboff;
	end

	% Data storage: each frame (one TR worth of data) is indexed by integers 'dabslice', 'dabecho', 'dabview'.
	% For now, only use the 'dabview' index. 
	dabslice = 0;
	dabecho  = 0;

	textra = 0;   % us
	freq = 0;     % Hz (RF offset frequency)
   waveform = 1; % waveform index (number; usually 1) -- each .mod file can have multiple different waveforms (of equal length)
	d = [d; ic ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview dabmode phi irfphase irfphase textra freq waveform];
end

% write matrix d to scanloop.txt
loop2txt(d);

return;


function g = g2gcm(g)
% convert gradient from Hz/cm to Gauss/cm

gamma = 4.2576e7;   % Hz/T
g = 1e2 * g / gamma;   % Gauss/cm

return;




%% OLD CODE %%

nt = size(d,1);              % number of startseq() calls   
NL = 13;   % toppe3.e 
maxslice = max(d(:,7));
maxecho = max(d(:,8));
maxview = max(d(:,9));

fname = 'scanloop.txt';
fid = fopen(fname, 'w', 'ieee-be');

fprintf(fid, 'nt\tmaxslice\tmaxecho\tmaxview\n');
fprintf(fid, '%d\t%d\t%d\t%d\n', nt, maxslice, maxecho, maxview);

fprintf(fid, 'Core ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview dabon phi rfphase recphase \n');
for ii = 1:nt
	for jj = 1:(NL-1)
		fprintf(fid, '%d\t', d(ii,jj));
	end
	fprintf(fid, '%d\n', d(ii,NL));
end

fclose(fid);

fprintf(1, 'done.\n');

return;


