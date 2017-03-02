function mat2mod(rho,theta,gx,gy,gz,nomflip,ofname,desc,hdrfloats,hdrints,addrframp,gsliceselect)
% function mat2mod(rho,theta,gx,gy,gz,nomflip,ofname,desc,[hdrfloats,hdrints,addrframp,gsliceselect])
% Write waveforms to .mod file, use with toppe.e psd on GE scanners.
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
% rho      -- Gauss
% theta    -- [-pi, pi]
% gx/gy/gz -- Gauss/cm
% ofname   -- output filename (e.g., mypulse.mod)
% desc     -- ASCII descriptor
% hdrfloats -- additional floats to put in header (max 12) (starting at paramsfloat[20])
% hdrints  -- [npre ndur [npre2 ndur2 npre3 ndur3 ...]] -- number of samples before start of RF/DAQ, 
%             and number of samples in DAQ (for one or more echoes)
% gsliceselect -- [G/cm] needed for RF frequency offset (coronal/sagittal scans)
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
% (c) 2017 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

if ~isreal(rho)
	exit('rho must be real valued');
end

if ~exist('addrframp','var')
	addrframp = 0;
end
if ~exist('gsliceselect','var')
	gsliceselect = 1;   % dummy value
end

% set gradient directions; add small RF ramp 
[rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp);

% calculate rfpulse info
[paramsfloat] = sub_myrfstat(rho(:), nomflip, gsliceselect);
if exist('hdrfloats','var')
	paramsfloat(20:(19+length(hdrfloats))) = hdrfloats;  % populate header with custom floats 
end
if exist('hdrints','var')
	paramsint16 = hdrints; 
else
	paramsint16 = [0 length(gx)];
end

% write to .mod file
desc = sprintf('Filename: %s\n%s', ofname, desc);
sub_writemod(ofname,desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat);

return;




function [rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp); %,tipupdelx,tipupdely)

%b1 = b1*1e4;     % convert to Gauss

% make smooth RF ramp to avoid RF amp error
ramp = rho(1)*[linspace(0,1,5)]';
ramp = [-ramp/2; flipud(-ramp/2); ramp];
if ~addrframp
	ramp = [0; 0];   % to avoid non-zero gradients at beginning/end which causes problems with readmod
end

rho = [ramp; rho];
theta = [0*ramp; theta];
gx = [0*ramp; gx];
gy = [0*ramp; gy];
gz = [0*ramp; gz];

% make length even
rho = makeevenlength(rho);   % units: Tesla
theta = makeevenlength(theta);   % units: Tesla
gx = makeevenlength(gx);   % units: G/cm
gy = makeevenlength(gy);
gz = makeevenlength(gz);

% fix gradient directions
%gx = -gx;
%gy = -gy;
%b1 = i*b1;

% Conjugate to compensate for modulator phase conjugation. Confirmed Feb 2011.
% headcoil: +phase on theta channel produced -phase in image -> also conjugate for headcoil
%b1 = conj(b1);

return;


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

if max(abs(b1))>0   % non-zero RF pulse
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
else % RF pulse is zero. This avoids division by zero.
	paramsfloat = [pw 0 0 0 0 ...
               0 1 0 0 ...
               0 nom_pw    nom_bw   1 0];
               0 

end
               
return;





function sub_writemod(fname,desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat)
%
%   desc           ASCII description 
%   rho            (Gauss) 
%                  size(rho) = N x Ncoils, where N = # waveform samples, 
%                  and Ncoils = number of transmit channels
%   theta          (radians). size(theta) = size(rho)
%   gx,gy,gz       (Gauss/cm). vector of length N
%   gmax           (Gauss/cm) maximum gradient amplitude supported by scanner  [4.0 G/cm]
%   paramsint16    Vector containing various int parameters, e.g., [npre ndur [additional optional values]]. Max length = 32.
%   paramsfloat    Vector containing various RF pulse parameters needed for SAR and B1 scaling
%                  calculations -- just a placeholder for now. Max length = 32.
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu

% max length of params* vectors
nparamsint16 = 32;
nparamsfloat = 32;

if max(abs(rho)) > 0.125
	max_b1 = max(abs(rho));
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

% are waveforms the right size?
rhosize   = size(rho);
thetasize = size(theta);
L = [max(rhosize) max(thetasize) numel(gx) numel(gy) numel(gz)];
if any(L-max(rhosize)) 
  error('mat2mod:unequallength', 'all waveforms must have the same length');
end
if any(rhosize ~= thetasize) 
  error('mat2mod:unequallength', 'rf and theta waveforms must have the same size and dimensions');
end

% convert waveforms to column-major form
[nr,nc] = size(rho);    if nc > nr;   rho   = rho';     end;
[nr,nc] = size(theta);  if nc > nr;   theta = theta';   end;
[nr,nc] = size(gx);     if nc > nr;   gx    = gx';      end;
[nr,nc] = size(gy);     if nc > nr;   gy    = gy';      end;
[nr,nc] = size(gz);     if nc > nr;   gz    = gz';      end;

% remember to make sure waveforms have even length

% convert params* to row vectors, and pad to max length
paramsint16  = reshape(paramsint16,   1, numel(paramsint16));
paramsfloat  = reshape(paramsfloat, 1, numel(paramsfloat));
paramsint16  = [paramsint16   zeros(1, nparamsint16-numel(paramsint16))];
paramsfloat  = [paramsfloat zeros(1, nparamsfloat-numel(paramsfloat))];

[res,ncoils] = size(rho);


%
% write to file
%
fid = fopen(fname, 'w', 'ieee-be');

% write header
globaldesc = sprintf('RF waveform file for ssfpbanding project.\n');  
globaldesc = sprintf('%sCreated by %s.m on %s.\n', globaldesc, mfilename('fullpath'), datestr(now));  
fs = dbstack;
if numel(fs)>2
	globaldesc = sprintf('%scalled by (dbstack(3)): %s\n', globaldesc, fs(3).file);
end
globaldesc = sprintf('%snchannels = %d, res = %d\n', globaldesc, ncoils, res);  
desc = sprintf('%s%s\n', globaldesc, desc);
fwrite(fid, numel(desc), 'int16');      % number of characters in ASCII description
fwrite(fid, desc, 'uchar');

fwrite(fid, ncoils, 'int16');          % shorts must be written in binary -- otherwise it won't work on scanner 
fwrite(fid, res,    'int16');
fprintf(fid, 'b1max:  %f\n', max_b1);   % floats are OK in ASCII on scanner 
fprintf(fid, 'gmax:   %f\n', gmax);

fwrite(fid, nparamsint16, 'int16');
fwrite(fid, paramsint16,  'int16');
fwrite(fid, nparamsfloat, 'int16');
for n = 1:nparamsfloat
	fprintf(fid, '%f\n', paramsfloat(n)); 
end

% write binary waveforms (*even* short integers -- the psd sets the EOS bit, so don't have to worry about it here)
maxiamp = 2^15;                           % max instruction amplitude (max value of signed short)
maxiamprf = 2^15-2;                       % RF amp is flipped if setting to 2^15 (as observed on scope), so stay below this
rho   = 2*round(rho/max_b1*maxiamprf/2);
theta = 2*round(theta/pi*maxiamp/2);
gx    = 2*round(gx/gmax*maxiamp/2);
gy    = 2*round(gy/gmax*maxiamp/2);
gz    = 2*round(gz/gmax*maxiamp/2);
fwrite(fid, rho(:),   'int16');
fwrite(fid, theta(:), 'int16');
fwrite(fid, gx(:),    'int16');
fwrite(fid, gy(:),    'int16');
fwrite(fid, gz(:),    'int16');

fclose(fid);

return;
