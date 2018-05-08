function [rho, th, gx, gy, gz] = scansim(nstart, nstop, d)
% function [rho, th, gx, gy, gz] = scansim(nstart, nstop, [d])
% display pulse sequence, as specified in cores.txt, scanloop.txt, and timing.txt
%
% $Id: scansim.m,v 1.2 2017/08/11 20:03:14 jfnielse Exp $

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
fid = fopen('modules.txt', 'r', 'ieee-be');
s = fgets(fid);  % skip line
ncores = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
for ic = 1:ncores
	cores{ic}.fname = fscanf(fid, '%s ', 1);
	cores{ic}.dur = fscanf(fid, '%d ', 1);
	cores{ic}.hasRF = fscanf(fid, '%d ', 1);
	cores{ic}.hasDAQ = fscanf(fid, '%d\n', 1);
	[desc,cores{ic}.rho,cores{ic}.th,cores{ic}.gx,cores{ic}.gy,cores{ic}.gz,cores{ic}.paramsint16,cores{ic}.paramsfloat] ...
		= readwav(cores{ic}.fname,false);
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

	rho = [rho; zeros(round((start_core+coredel)/dt),1); ia_rf/max_pg_iamp*cores{ic}.rho; zeros(round((timetrwait+timessi)/4),1);         zeros(round(textra/4),1)];
	th  = [th;  zeros(round((start_core+coredel)/dt),1); ia_th/max_pg_iamp*cores{ic}.th;  zeros(round((timetrwait+timessi)/4),1);         zeros(round(textra/4),1)];
	gx  = [gx;  zeros(round((start_core)/dt),1);         ia_gx/max_pg_iamp*cores{ic}.gx;  zeros(round((timetrwait+timessi+coredel)/4),1); zeros(round(textra/4),1)];
	gy  = [gy;  zeros(round((start_core)/dt),1);         ia_gy/max_pg_iamp*cores{ic}.gy;  zeros(round((timetrwait+timessi+coredel)/4),1); zeros(round(textra/4),1)];
	gz  = [gz;  zeros(round((start_core)/dt),1);         ia_gz/max_pg_iamp*cores{ic}.gz;  zeros(round((timetrwait+timessi+coredel)/4),1); zeros(round(textra/4),1)];

	%fprintf(1, 'it %d: t = %.3f ms\n', it, numel(rho)*dt*1e-3);
end

% plot
T = (1:numel(rho))*dt/1000; % msec
gmax = 5;  % Gauss/cm
subplot(511); plot(T, rho); ylabel('rho');
subplot(512); plot(T, th);  ylabel('theta');
subplot(513); plot(T, gx);  ylabel('gx'); axis([T(1) T(end) -gmax gmax]);
subplot(514); plot(T, gy);  ylabel('gy'); axis([T(1) T(end) -gmax gmax]);
subplot(515); plot(T, gz);  ylabel('gz'); axis([T(1) T(end) -gmax gmax]);
xlabel('msec');

return;


% EOF
