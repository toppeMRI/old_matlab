function d = writeloop
% Toy example: cycle through 12 different readout waveforms.
%
% INPUTS:
%   freq*: frequency offsets (Hz) of 90 and 180 pulses. Get this value from makepencil.m (see main.m).
%
% $Id: writeloop_pencil.m,v 1.26 2016/12/07 17:52:49 jfnielse Exp $

fprintf(1,'\twriting scanloop.txt...');

% acquisition parameters
nl = 2;           % Number of readout waveform in-plane rotations to cycle through
nave = 1;         % number of averages.
nwaveforms = 12;     % number of GIRF calibration waveforms

% define some constants
daboff = 0;
dabon  = 1;                  % for nex=1
dabadd = 2;                  % for nex > 1
max_pg_iamp = 2^15-2;

freq90 = 0;    % in an actual GIRF calibration, these would be non-zero
freq180 = 0;

ia_rf = max_pg_iamp;    % full scale (90 or 180)
ia_th = max_pg_iamp;

% play DAB reset core, then loop through acquisitions
d = [];
for ifr = 1:nave 
	dabslice = 1;
	if ifr==1
		dabmode = dabon;
	else
		dabmode = dabadd;
	end
	for ileaf = 1:nl
		dabview = ileaf;                               % skip baseline view as always  
		phi = round(2*pi*(ileaf-1)/nl/pi*max_pg_iamp);
		for iecho = 1:nwaveforms
			dabecho = iecho-1;
			irfphase  = 0;
			irecphase = 0;

			% 90 degree excitation
			core = 1; 
			waveform = 1;
			textra = 0;
			d = [d; core ia_rf ia_th max_pg_iamp 0*max_pg_iamp 0*max_pg_iamp 0 0 0 daboff 0 irfphase irecphase textra round(freq90) waveform]; 

			% 180 degree spin-echo pulse
			core = 2; 
			waveform = 1;
			textra = 0;
			d = [d; core ia_rf ia_th 0*max_pg_iamp max_pg_iamp 0*max_pg_iamp 0 0 0 daboff 0 irfphase irecphase textra round(freq180) waveform]; 

			% readout
			core = 3; 
			waveform = iecho;
			textra = 0;
			ia_gx = max_pg_iamp;
			ia_gy = max_pg_iamp;
			d = [d; core 0 0 ia_gx ia_gy 0 dabslice dabecho dabview dabmode phi irfphase irecphase textra 0 waveform];
		end
	end
end

% write loop file
loop2txt(d);

return;
