function dur = getScanTime(varargin)
% function dur = getScanTime(varargin)
%
% Estimate total scan time.
%
% Input options:
%   'scanloopfile'     default: 'scanloop.txt'
%   'modulesfile'      default: 'modules.txt'
%   'timingfile'       default: 'timing.txt'
%   In addition, the .mod files listed in 'modulesfile' must be present in the current (working) folder.
%
% Output:
%    dur          sec
%

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
% (c) 2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: getScanTime.m,v 1.6 2018/10/08 14:23:20 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2/getScanTime.m,v $

%% parse inputs
% Default values 
arg.scanloopfile = 'scanloop.txt';
arg.modulesfile = 'modules.txt';
arg.timingfile = 'timing.txt';

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

%% read files
% read scanloop
looparr = tryread(@readloop, arg.scanloopfile);

% read timing parameters (CVs)
TPARAMS = tryread(@readTimingFile, arg.timingfile);

% read module content
mods = tryread(@readModules, arg.modulesfile);

%% loop through scan, and tally scan duration
dt = 4e-6;    % duration of one gradient/rf sample (sec)
dur = 0;
fprintf('Looping through scan and tallying scan duration...');
for ii = 1:size(looparr,1)
	if ~mod(ii,10000)
		fprintf('.');
	end
	rho = dispseq(ii, ii, 'looparr', looparr, 'tparams', TPARAMS, 'mods', mods, 'dodisplay', false);
	dur = dur + size(rho,1)*dt;
end
fprintf(' done\n');
dur = round(dur);
fprintf('Total scan time: %dm %ds\n', floor(dur/60), dur - 60*floor(dur/60) );

