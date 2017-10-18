function dat = loadpfile(pfile)
% function dat = loadpfile(pfile)
%
% Load all data from Pfile, EXCEPT dabslice=0 slot (which can contain corrupt data).

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

% read Pfile header
fid = fopen(pfile,'r','l');
rdb_hdr = read_rdb_hdr(fid,24);
fclose(fid);

ndat    = rdb_hdr.frame_size;
nslices = rdb_hdr.nslices;
nechoes = rdb_hdr.nechoes;
nviews  = rdb_hdr.nframes;
ncoils  = rdb_hdr.dab(2)-rdb_hdr.dab(1)+1;

%dat = zeros([ndat ncoils nslices nechoes nviews]);
fid = fopen(pfile,'r','l');
for slice = 2:nslices   % skip first slice (sometimes contains corrupted data)
	for echo = 1:nechoes
		for view = 1:nviews
			[dattmp pfilesize] = loaddat_ge(fid,rdb_hdr,slice-1,echo-1,view);     % [ndat ncoils]. Skip baseline (0) view.
			dat(:,:,slice,echo,view) = dattmp; %(:,1:16:ncoils);                                     
		end
	end
end
fclose(fid);
fprintf(1,'Expected pfilesize = %d\n', pfilesize);
return;

% Average
%dat = mean(dat,5);

% Subtract reference scan phase.
%dat = dat(:,:,:,1:4)./exp(1i*angle(dat(:,:,:,5:8)));          % [ndat ncoils nleafs 4 nframes] Subtract reference frame.


% Combine data phase from all coils, weight coils by intensity squared.
datw = 0;
for ic = 1:ncoils
	d1 = squeeze(dat(:,ic,:,:,:));               
	datw = datw + abs(d1).^2.*exp(1i*angle(d1));             % [ndat nleafs 4]
end
th = angle(datw);                                           % [ndat nleafs 4]

return;
