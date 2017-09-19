function playseq(loopfile,nModPerTR,nTRskip,tpause)
% function playseq(loopfile,nModPerTR,nTRskip)
%
% Loop through a scanloop.txt file ("play" sequence)
%
% INPUTS:
%    loopfile:   'scanloop.txt'
%    nModPerTR:  number of modules per TR
%    nTRskip:    display only every nTRskip TRs (for speeding up loop)
%    tpause:     delay before next TR (sec)
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
% (c) 2016 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: playseq.m,v 1.1 2017/02/15 18:24:09 jfnielse Exp $

if ~exist('tpause','var')
	tpause = 0.02;   % sec
end

d = readloop(loopfile);
nt = size(d,1);

for ii = 1:((1+nTRskip)*nModPerTR):nt
	dispseq(ii,ii+nModPerTR-1,d);
	title(num2str(ii));
	pause(tpause);     % to allow display to refresh
end
