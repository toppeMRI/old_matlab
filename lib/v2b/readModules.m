function mods = readModules(modulesfile)
% function mods = readModules(modulesfile)
%
% $Id: readModules.m,v 1.2 2018/10/09 15:22:34 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2/readModules.m,v $

if ~exist('modulesfile', 'var')
	modulesfile = 'modules.txt';
end

% get waveforms
fid = fopen(modulesfile, 'r', 'ieee-be');
s = fgets(fid);  % skip line
ncores = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
for ic = 1:ncores
	cores{ic}.fname = fscanf(fid, '%s ', 1);
	cores{ic}.dur = fscanf(fid, '%d ', 1);
	cores{ic}.hasRF = fscanf(fid, '%d ', 1);
	cores{ic}.hasDAQ = fscanf(fid, '%d\n', 1);
	[desc,cores{ic}.rf,cores{ic}.gx,cores{ic}.gy,cores{ic}.gz,cores{ic}.paramsint16,cores{ic}.paramsfloat] ...
		= readmod(cores{ic}.fname,false);
	cores{ic}.wavdur = numel(cores{ic}.gx(:,1))*4;   % waveform duration [us]
end
fclose(fid);

mods = cores;
