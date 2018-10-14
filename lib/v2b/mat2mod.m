function mat2mod(varargin)
% function mat2mod(varargin)
%
% Write waveforms to .mod file, use with toppev2 psd on GE scanners.
% Assumes raster time (sample duration) of 4e-6 sec for all waveforms.
%
% Examples:
% >> mat2mod('rf', rho.*exp(1i*theta), 'gz', gzwaveform, 'nomflip', 30);
% >> mat2mod('gz', gzwaveform, 'desc', 'my spoiler gradient');
% >> lims = {'MaxGrad', 50,  'GradUnit', 'mT/m', ...
%            'MaxSlew', 200, 'SlewUnit', 'T/m/s'};
% >> mat2mod('rf, myrf, 'gx', gzwav, 'system', lims);
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
% Input options:
%   rf            Complex RF waveform, [ndat nrfpulses]
%   gx            [ndat ngxpulses]
%   gy            [ndat ngypulses]
%   gz            [ndat ngzpulses]
%   rfUnit        'Gauss' (default) or 'mT'
%   gradUnit      'Gauss/cm' (default) or 'mT/m'
%   ofname        Output filename.
%   desc          Text string (ASCII) descriptor.
%   nomflip       Excitation flip angle (degrees). Default: 90.
%   hdrfloats     Additional floats to put in header (max 12)
%   hdrints       Additional ints to put in header (max 30)
%   system        Cell array specifying hardware system limits (per manufacturer specification).
%                 Default:
%                 system = {'MaxGrad', 50,  'GradUnit', 'mT/m', ...
%                           'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
%                           'MaxRf',   20,  'RfUnit',   'mT'};  
%
% $Id: mat2mod.m,v 1.4 2018/10/14 14:54:16 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2b/mat2mod.m,v $

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
% (c) 2017-18 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

%% parse inputs
% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];
arg.rfUnit    = 'Gauss';
arg.gradUnit  = 'Gauss/cm';
arg.ofname    = 'out.mod'; 
arg.desc      = 'TOPPE module';
arg.nomflip   = 90;
arg.hdrfloats = [];
arg.hdrints   = [];
arg.system    = [];

arg = vararg_pair(arg, varargin);

% parse system hardware limits
% First, define defaults.
system.MaxGrad  = 50;
system.GradUnit = 'mT/m';
system.MaxSlew  = 200;
system.SlewUnit = 'T/m/s';
system.MaxRf    = 20;
system.RfUnit   = 'mT';
if ~isempty(arg.system)
	% Substitute specified system values as appropriate
	system = vararg_pair(system, arg.system);
end

%% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
	wavtype = fields{ii} ;   % 'rf', 'gx', 'gy', or 'gz'
	cmd = sprintf('%s = %s;', wavtype, sprintf('arg.%s', wavtype));
	eval(cmd);
end

%% Convert to Gauss and Gauss/cm
if strcmp(arg.rfUnit, 'mT')
	rf = rf/100;   % Gauss
end
if strcmp(arg.gradUnit, 'mT/m')
	gx = gx/10;    % Gauss/cm
	gy = gy/10;
	gz = gz/10;
end


%% Force all waveform arrays to have the same dimensions (required by toppev2.e)
ndat    = max( [size(rf,1) size(gx,1) size(gy,1) size(gz,1)] );
npulses = max( [size(rf,2) size(gx,2) size(gy,2) size(gz,2)] );
if ~ndat 
	error('At least one waveform must be specified');
end

for ii = 1:length(fields)
	wavtype = fields{ii} ;                    % 'rf', 'gx', 'gy', or 'gz'
	wav = eval(fields{ii}) ;                  % [ndat npulses]

	if wavtype == 'rf' & isempty(wav)
		wav = 0.01*ones(ndat, npulses);        % Must have non-zero RF waveform to make toppev2.e happy (even if it's not actually played out)
	end

	% enforce equal number of rows and columns
	[nrows n2] = size(wav);

	if (nrows ~=0 & nrows < ndat) 
		warning('Padding %s with zero rows', wavtype);
	end
	if (n2 ~= 0 & n2 < npulses) 
		warning('Padding %s with zero columns', wavtype);
	end

	wav = [wav; zeros(ndat-nrows,n2)];
	wav = [wav  zeros(ndat,npulses-n2)];

	wav = [zeros(2,npulses); wav];   % to avoid non-zero gradients at beginning which causes problems with readmod

	% make length divisible by 4 (EPIC seems to behave best this way)
	wav = makeGElength(wav);

	% copy to corresponding wav type (rf, gx, gy, or gz)
	cmd = sprintf('%s = %s;', wavtype, 'wav') ;
	eval (cmd);
end

%% Check against system hardware limits
if strcmp(system.RfUnit, 'mT')
	system.MaxRf = system.MaxRf/100;      % Gauss
end
if strcmp(system.GradUnit, 'mT/m')
	system.MaxGrad = system.MaxGrad/10;   % Gauss/cm
end
if strcmp(system.SlewUnit, 'T/m/s')
	system.MaxSlew = system.MaxSlew/10;   % Gauss/cm/msec
end

systemViolation = false;

maxG = max([abs(gx(:)) abs(gy(:)) abs(gz(:))]);
grads = 'xyz';
for ii = 1:3
	if maxG(ii) > system.MaxGrad
		fprintf('Error: %s gradient amplitude exceeds system limit (%.1f%%)\n', grads(ii), maxG(ii)/system.MaxGrad*100);
		systemViolation = true;
	end
end

rasterTime = 4e-3;
gslew = diff([gx gy gz]/rasterTime,1);
maxslew = max(gslew);
grads = 'xyz';
for ii = 1:3
	if maxslew(ii) > system.MaxSlew
		fprintf('Error: %s gradient slew rate exceeds system limit (%.1f%%)\n', grads(ii), maxslew(ii)/system.MaxSlew*100);
		systemViolation = true;
	end
end

maxRf = max(rf);
if maxRf > system.MaxRf
	fprintf('Error: rf amplitude exceeds system limit (%.1f%%)\n', maxRf/system.MaxRf*100);
	systemViolation = true;
end

if systemViolation
	error('Errors found -- exiting');
end
	

% Fixes to avoid idiosyncratic issues on scanner
%[rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp);

%% Optional header arrays
gsliceselect = 1;  % legacy dummy value, ignore
[paramsfloat] = sub_myrfstat(abs(rf(:,1,1)), arg.nomflip, gsliceselect);
if ~isempty(arg.hdrfloats)
	paramsfloat(20:(19+length(arg.hdrfloats))) = arg.hdrfloats;  % populate header with custom floats 
end
paramsint16 = [0 size(rf,1)];
if ~isempty(arg.hdrints)
	paramsint16(3:(2+length(arg.hdrints))) = arg.hdrints;  % populate header with custom ints
end

%% Write to .mod file
arg.desc = sprintf('Filename: %s\n%s', arg.ofname, arg.desc);
sub_writemod(arg.ofname, arg.desc, rf, gx, gy, gz, paramsint16, paramsfloat);

return;



%% function paramsfloat = sub_myrfstat(b1, nom_fa, g);
function paramsfloat = sub_myrfstat(b1, nom_fa, g);
% Calculate RF parameters needed for RFPULSE struct in .e file.
% Needed for B1 scaling, SAR calculations, and enforcing duty cycle limits.
% See also mat2signa_krishna.m
%
% b1         real 1D vector containing B1 amplitude, size Nx1 [Gauss]
% nom_fa     nominal flip angle (degrees)

nom_bw = 2000;

dt = 4e-6;                        % use 4 us RF sample width
%gamma = 4.2575e3;                  % Hz/Gauss
tbwdummy = 2;

hardpulse = max(abs(b1)) * ones(length(b1),1);    % hard pulse of equal duration and amplitude

pw            = length(b1)*dt*1e3;                                       % ms
nom_pw        = length(b1)*dt*1e6;                                        % us

if max(abs(b1)) == 0   % non-zero RF pulse
	error('RF waveform cannot be zero');
end
abswidth      = sum(abs(b1)) / sum(abs(hardpulse));
effwidth      = sum(b1.^2)   / sum(hardpulse.^2);
% or equivalently:  effwidth = sum(b1.^2)/(max(abs(b1))^2)/length(b1)
area          = abs(sum(b1)) / abs(sum(hardpulse)); 
dtycyc        = length(find(abs(b1)>0.2236*max(abs(b1)))) / length(b1);
maxpw         = dtycyc;
num           = 1;
if max(abs(b1)) > 0.125
	max_b1        = max(abs(b1)) ;                                            % Gauss
else
	max_b1        = 0.125 ;                                                     % Gauss
end
max_int_b1_sq = max( cumsum(abs(b1).^2)*dt*1e3 );                              % Gauss^2 - ms
max_rms_b1    = sqrt(mean(abs(b1).^2));                                        % Gauss
nom_fa        = nom_fa;                                                   % degrees
%nom_bw        = tbwdummy / (dt * length(b1));                             % Hz
% max_int_b1    = abs(sum(b1))*dt*1000

% calculate equivalent number of standard pulses
stdpw = 1;                                                % duration of standard pulse (ms)
stdpulse = 0.117 * ones(round(stdpw/(dt*1e3)),1);
numstdpulses = num * effwidth * (pw/stdpw) * (max(abs(b1))/0.117)^2;

%pulse50 = 0.117/180*50 * ones(round(stdpw/(dt*1e3)),1);
%num50 = sum(abs(b1).^2)/sum(pulse50.^2);

paramsfloat = [pw     abswidth  effwidth area          dtycyc      ...
              maxpw  num       max_b1   max_int_b1_sq max_rms_b1  ...
              nom_fa nom_pw    nom_bw   g numstdpulses              ];

%else % RF pulse is zero. This avoids division by zero.
%	paramsfloat = [pw 0 0 0 0 ...
%               0 1 0 0 ...
%               0 nom_pw    nom_bw   1 0];
%               0 
%end
               
return;



%% function sub_writemod(fname,desc,rf,gx,gy,gz,paramsint16,paramsfloat)
function sub_writemod(fname,desc,rf,gx,gy,gz,paramsint16,paramsfloat)
%
%   desc           ASCII description 
%   rf             size(rho) = N x 1, where N = # waveform samples, 
%                  and Ncoils = number of transmit channels
%   gx,gy,gz       (Gauss/cm). vector of length N
%   paramsint16    Vector containing various int parameters, e.g., [npre ndur [additional optional values]]. Max length = 32.
%   paramsfloat    Vector containing various RF pulse parameters needed for SAR and B1 scaling
%                  calculations -- just a placeholder for now. Max length = 32.
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu

rho = abs(rf);
theta = angle(rf);

npulses = size(rf,2);

% max length of params* vectors
nparamsint16 = 32;
nparamsfloat = 32;

if max(abs(rho)) > 0.125
	max_b1 = max(abs(rho(:)));
else
	max_b1 = 0.125;      % waveform is scaled relative to this. Max b1 is 0.25G/0.125G for quadrature/body RF coils (according to John Pauly RF class notes)
end

gmax  = 5.0;                  % Gauss/cm
b1max = max(abs(rho(:)));     % Gauss

if (numel(paramsint16)>nparamsint16)
  error('mat2mod:nparamsint16',   'too many int16 parameters');
end
if (numel(paramsfloat)>nparamsfloat)
  error('mat2mod:nparamsfloat', 'too many float parameters');
end

% convert params* to row vectors, and pad to max length
paramsint16  = reshape(paramsint16,   1, numel(paramsint16));
paramsfloat  = reshape(paramsfloat, 1, numel(paramsfloat));
paramsint16  = [paramsint16   zeros(1, nparamsint16-numel(paramsint16))];
paramsfloat  = [paramsfloat zeros(1, nparamsfloat-numel(paramsfloat))];

[res,nrfpulses,ncoils] = size(rho);

fid = fopen(fname, 'w', 'ieee-be');

% write header
globaldesc = sprintf('RF waveform file for ssfpbanding project.\n');  
globaldesc = sprintf('%sCreated by %s.m on %s.\n', globaldesc, mfilename('fullpath'), datestr(now));  
fs = dbstack;
if numel(fs)>2
	globaldesc = sprintf('%scalled by (dbstack(3)): %s\n', globaldesc, fs(3).file);
end
globaldesc = sprintf('%sncoils = %d, res = %d\n', globaldesc, ncoils, res);  
desc = sprintf('%s%s\n', globaldesc, desc);
fwrite(fid, numel(desc), 'int16');      % number of characters in ASCII description
fwrite(fid, desc, 'uchar');

fwrite(fid, ncoils,  'int16');          % shorts must be written in binary -- otherwise it won't work on scanner 
fwrite(fid, res,     'int16');
fwrite(fid, npulses, 'int16');
fprintf(fid, 'b1max:  %f\n', max_b1);   % floats are OK in ASCII on scanner 
fprintf(fid, 'gmax:   %f\n', gmax);

fwrite(fid, nparamsint16, 'int16');
fwrite(fid, paramsint16,  'int16');
fwrite(fid, nparamsfloat, 'int16');
for n = 1:nparamsfloat
	fprintf(fid, '%f\n', paramsfloat(n)); 
end

% write binary waveforms (*even* short integers -- the toppe.e psd sets the EOS bit, so don't have to worry about it here)
max_pg_iamp = 2^15-2;              % RF amp is flipped if setting to 2^15 (as observed on scope), so subtract 2
rho   = 2*round(rho/max_b1*max_pg_iamp/2);
theta = 2*round(theta/pi*max_pg_iamp/2);
gx    = 2*round(gx/gmax*max_pg_iamp/2);
gy    = 2*round(gy/gmax*max_pg_iamp/2);
gz    = 2*round(gz/gmax*max_pg_iamp/2);

for ip = 1:npulses
	for ic = 1:ncoils
		fwrite(fid, rho(:,ip,ic), 'int16');
	end
	for ic = 1:ncoils
		fwrite(fid, theta(:,ip,ic),   'int16');
	end
	fwrite(fid, gx(:,ip),    'int16');
	fwrite(fid, gy(:,ip),    'int16');
	fwrite(fid, gz(:,ip),    'int16');
end

fclose(fid);

return;



%% function [rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp); %,tipupdelx,tipupdely)
function [rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp); %,tipupdelx,tipupdely)

% make smooth RF ramp to avoid RF amp error
ramp = rho(1)*[linspace(0,1,5)]';
ramp = [-ramp/2; flipud(-ramp/2); ramp];
if ~addrframp
	ramp = [0; 0];   % to avoid non-zero gradients at beginning/end which causes problems with readmod
end

[n npulses ncoils] = size(rho);   % ncoils is number of RF transmit coils 
for ip = 1:npulses
	for ic = 1:ncoils
		rho2(:,ip,ic) = [ramp; rho(:,ip,ic)];
		theta2(:,ip,ic) = [0*ramp; theta(:,ip,ic)];
	end
end
rho = rho2;
theta = theta2;
gx = [zeros(numel(ramp),npulses); gx];
gy = [zeros(numel(ramp),npulses); gy];
gz = [zeros(numel(ramp),npulses); gz];


return;
