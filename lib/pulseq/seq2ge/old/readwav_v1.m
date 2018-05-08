function [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav(fname,showinfo)
% function [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav(fname,showinfo)
%
% Read .wav file
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu
% $Id: readwav.m,v 1.1 2015/07/01 12:46:20 jfnielse Exp $

if ~exist('showinfo','var')
	showinfo = true;
end

fid = fopen(fname, 'r', 'ieee-be');

% read ASCII description
asciisize = fread(fid, 1,         'int16');
desc      = fread(fid, asciisize, 'uchar'); 
desc = char(desc');

% read rest of header
ncoils = fread(fid, 1, 'int16');
res    = fread(fid, 1, 'int16');
b1max  = fscanf(fid, 'b1max:  %f\n');
gmax   = fscanf(fid, 'gmax:   %f\n');

nparamsint16 = fread(fid, 1,            'int16');
paramsint16  = fread(fid, nparamsint16, 'int16');
nparamsfloat = fread(fid, 1,            'int16');
%paramsfloat  = fscanf(fid, '%f\n');                % this reads all float32 params at once
for n = 1:nparamsfloat
	paramsfloat(n)  = fscanf(fid, '%f\n', 1);
end

%rho = []; theta = []; gx = []; gy = []; gz = []; return;

if showinfo
	fprintf(1, '\n%s', desc);
	fprintf(1, 'number of coils/channels:  %d\n', ncoils);
	fprintf(1, 'number points in waveform: %d\n', res);
	fprintf(1, 'data offset (bytes):       %d\n', ftell(fid));
	fprintf(1, '\n');
end

% read waveforms (short)
rho   = fread(fid, res*ncoils, 'int16');
theta = fread(fid, res*ncoils, 'int16');
gx    = fread(fid, res,        'int16');
gy    = fread(fid, res,        'int16');
gz    = fread(fid, res,        'int16');

% convert back to physical units
maxiamp = 2^15-2;                           % max instruction amplitude (max value of signed short)
rho   = rho*b1max/maxiamp;     % Gauss
theta = theta*pi/maxiamp;      % radians
gx    = gx*gmax/maxiamp;       % Gauss/cm
gy    = gy*gmax/maxiamp;
gz    = gz*gmax/maxiamp;

% reshape output
rho   = reshape(rho,   res, ncoils);
theta = reshape(theta, res, ncoils);
gx    = reshape(gx,    res, 1);
gy    = reshape(gy,    res, 1);
gz    = reshape(gz,    res, 1);
paramsint16 = paramsint16';
paramsfloat = paramsfloat';

fclose(fid);

% EOF
