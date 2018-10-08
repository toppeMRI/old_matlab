function playseq(loopfile,nModPerTR,varargin)
% function playseq(loopfile,nModPerTR,varargin)
%
% Loop through a scanloop.txt file ("play" sequence)
%
% INPUTS:
%    loopfile:   'scanloop.txt'
%    nModPerTR:  number of modules per TR
% Options:
%    nTRskip:      display only every nTRskip TRs (for speeding up loop) (default: 0)
%    tpause:       delay before displaying next TR (sec) (default: 0)
%    modulesfile   Default: 'modules.txt'
%    timingfile    Default: 'timing.txt'
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
% (c) 2016-18 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: playseq.m,v 1.2 2018/10/08 14:14:01 jfnielse Exp $

if nargin < 2
	error('nModPerTR must be specified');
end
if nargin < 1
	error('Scan loop file and nModPerTR must be specified');
end

%% parse inputs
% Default values 
arg.nTRskip = 0;
arg.tpause = 0.02; 
arg.mods = [];
arg.modulesfile = 'modules.txt';
arg.tparams = [];
arg.timingfile = 'timing.txt';

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

%% read scan files as needed

% read scanloop
looparr = tryread(@readloop, loopfile);

% get timing CVs
if isempty(arg.tparams)
	TPARAMS = tryread(@readTimingFile, arg.timingfile);
else
	TPARAMS = arg.tparams;
end
[start_core myrfdel daqdel timetrwait timessi] = deal(TPARAMS(1), TPARAMS(2), TPARAMS(3), TPARAMS(4), TPARAMS(5));

% read module waveforms
if isempty(arg.mods)
	cores = tryread(@readModules, arg.modulesfile);
else
	cores = arg.mods;
end

%% display
for ii = 1 : ((1+arg.nTRskip)*nModPerTR) : size(looparr,1)
	dispseq(ii, ii+nModPerTR-1, 'looparr', looparr);
	subplot(511); title(num2str(ii));
	pause(arg.tpause);     % to allow display to refresh
end
