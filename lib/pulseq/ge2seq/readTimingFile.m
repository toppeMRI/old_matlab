function TPARAMS = readTimingFile(timingfile)
% function TPARAMS = readTimingFile([timingfile])
%
% Inputs:
%   timingfile	        Default: 'timing.txt'
%
% Outputs:
%   TPARAMS     [start_core myrfdel daqdel timetrwait timessi] (1x5 array)

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
% (c) 2016-2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: readTimingFile.m,v 1.1 2018/10/09 22:17:56 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/pulseq/pulseq-master/matlab/ge2seq/readTimingFile.m,v $

if ~exist('timingfile', 'var')
	timingfile = 'timing.txt';
end

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
	
