function d = writeloop
%
% Stack-of-spirals whole-brain 3D PRESTO fMRI sequence, 3.33mm isotropic.
% Acquire a few fully sampled frames at beginning (for reference and B0 mapping), 
% then undersample (R~6).
% Undersampling pattern is the same for each frame, which enables ghost removal by 
% temporal filtering (Nielsen and Noll, MRM 2016).
%
% $Id: writeloop.m,v 1.37 2017/02/14 21:01:58 jfnielse Exp $

% acquisition parameters
% 3.33 mm isotropic resolution
%rhrecon = 810;  % Stack-of-spirals recon script is /usr/g/bin/recon810
nleafs = 3;     % Number of spiral rotations (leafs).
nz = 54;        % Number of fully sampled z phase-encodes. Excitation is 14 cm, readout FOV in z is 54*3.33mm = 180mm
nz_samp = 30;   % Sample this many kz points per time-frame
                % NB! To make ghost phase alternate b/w 0 and 180 (relative to main object) ever other frame, 
                % we want mod(nleafs*nz_samp*rf_spoil_seed,360) = 180. 
trvol = 3;      % (approximate) time for fully sampled image volume (sec)
dur = 4*60;     % total duration of rs-fMRI scan (sec)

%rf_spoil_seed = 150;    % for nleafs=3, nz=60, and undersampling Rxy*rz=6, we have 30 TRs per time-frame.
rf_spoil_seed = 150;    % for 30 kz platters per frame and rf_spoil_seed=150, we have mod(nz_samp*rf_spoil_seed)=180.

% undersampling factors
Rxy = 3;
Rz = nz/nz_samp;

nt = 2*round(Rxy*Rz*dur/trvol/2);      % number of (undersampled) time-frames

% define some constants
dabon = 1;
daboff = 0;
max_pg_iamp = 2^15-2;

% set flip angle
nomflipd = 90;   % nominal (design) tipdown angle
ia_rf = 2*round([10]/nomflipd*max_pg_iamp/2);
ia_th = max_pg_iamp;

% fully sampled kz sampling pattern
for ii = 1:nz
	kzFull(ii) = ((ii-1+0.5)-nz/2)/(nz/2);    % scaling is (-1 1)
end

if 0
	% Variable-density (non-Cartesian) kz undersampling
	kz = vardenskz(nz,Rz,3.3);    % Fully sampled center with quadratically increasing FOV outside center. Last argument is FOV(kz=0)/FOV(kzmax). 
	a_gz_max = abs((0.5-nz/2)/(nz/2));
	kzU = kz*a_gz_max;     % scaling is (-1 1)
else
	% Cartesian variable-density kz undersampling
	load zInd;
	kzU = kzFull(logical(zInd));
end

% get spiral leaf
%[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav('readout.mod');
%[kinfo.kx kinfo.ky] = g2k([gx(:) gy(:)],nleafs);

% loop through acquisitions
d = []; 
dabecho = 0;
irfphase = 0;
rfphase = 0;
rf_spoil_seed_cnt = 0;

textra = 0;

ndisdaq = 2;
nref = 4*nleafs;     % (number of fully sampled frames acquired at beginning) * (nleafs)

% initialize d array (faster this way)
NL = 15;             % Number of columns in scanloop.txt (toppe8*)
d = zeros((ndisdaq+nref+nt)*nz*2,NL);   % 2 startseq() calls per TR. Will be cropped to right size after loop.

% Loop through acquisitions.
% At each frame, kxy is undersampled by factor 3.
% To get a fully sampled reference set, combine three consecutive reference frames.
icnt = 1;      % counts startseq() calls
for iframe = (-ndisdaq+1):(nref+nt)

	% Set PRESTO gradient amplitude.
	% Step through a few different amplitudes.
	nsteps = 4;
	nFramesPerStep = round(nt/nsteps);
	if iframe < (nref+1)    % disdaqs or reference frames
		a_presto = 1;
	else
		a_presto =  1 - (floor(( (iframe-nref)-1)/nFramesPerStep))*1/nsteps;
	end
	ia_presto = 2*round(a_presto*max_pg_iamp/2);
	
	% disdaqs
	if iframe < 1 
		dabmode = daboff;
		dabview  = 0;                    % this doesn't seem necessary but seems good practice
	else
		dabmode = dabon;
		dabview = iframe;
	end

	% Set kz sampling pattern for this frame.
	if iframe < (nref+1)
		kz = kzFull;                % numel(kz) = nz;
	else
		kz = kzU;                   % numel(kzU) = nz_samp;
	end

	for iz = 1:numel(kz)

		dabslice = iz;      % skip dabslice = 0

		%ia_gz = 2*round( max_pg_iamp*(((iz-1+0.5)-nz/2)/(nz/2)) /2);
		ia_gz = 2*round(kz(iz)*max_pg_iamp/2);

		if iframe < (nref+1)
			ileaf = mod(mod(iframe,nleafs) + iz,nleafs) + 1;   % rotate leaf every frame and every kz platter
		else
			ileaf = mod(iz,nleafs) + 1;                        % rotate leaf every kz platter (i.e., same undersampling pattern for every frame)
		end

		dabecho = 0 ;%ileaf-1;

		%phi = round(2*pi*(ileaf-1)/nleafs/pi*max_pg_iamp);
		phi = 2*pi*(ileaf-1)/nleafs;                            % [0 2pi]
		if phi > pi
			phi = phi - 2*pi;                                    % wrapped to [-pi pi]
		end    
		iphi = 2*round(phi/pi*max_pg_iamp/2);

		% transmit and receive phase
		irecphase = irfphase;   % phase from previous TR (PRESTO)
		rfphase  = rfphase + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
		rfphasetmp = atan2(sin(rfphase), cos(rfphase));       % wrap phase to (-pi,pi) range
		irfphase = 2*round(rfphasetmp/pi*max_pg_iamp/2);      % even short int 

		% slice-select (tipdown) pulse, and PRESTO gradients
		core = 1; 
		d(icnt,:) = [core ia_rf ia_th ia_presto ia_presto max_pg_iamp 0 0 0 daboff 0 irfphase irecphase textra 0]; 
		icnt = icnt+1;

		% spiral-in readout (balanced)
		core = 2; 
		d(icnt,:) = [core 0 0 max_pg_iamp max_pg_iamp ia_gz dabslice dabecho dabview dabmode iphi irfphase irecphase textra 0]; 
		%kinfo.kind(icnt,:) = [ileaf ia_gz];
		icnt = icnt+1;

	end
end

%disp('saving kinfo')
%kinfo.Rxy = Rxy;
%kinfo.Rz = Rz;
%kinfo.nref = nref;
%kinfo.ndisdaq = ndisdaq;
%save -v7.3 kinfo kinfo

d = d(1:(icnt-1),:);  % in case d is too big

% write loop file
loop2txt(d);

return;


function kzindU = sub_samp(nz,Rz,r)

kzindU = 1:2:nz;

return;



