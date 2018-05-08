function d = readloop(fname)
% function d = readloop(fname)
% read scanloop.txt (toppe.e scan loop definition, used with cores.txt).
% $Id: readloop.m,v 1.2 2017/02/07 16:31:19 jfnielse Exp $

d = []; 
NL = 14;   % toppe5
NL = 11;   % toppe2
NL = 10;   % toppe.e and toppe_so.e
NL = 13;   % toppe3, toppe4

fid = fopen(fname, 'r', 'ieee-be');

s = fgets(fid);  % read a whole line

nt       = fscanf(fid, '%d\t', 1);
maxslice = fscanf(fid, '%d\t', 1);
maxecho  = fscanf(fid, '%d\t', 1);
maxview  = fscanf(fid, '%d\t', 1);
rhrecon = fscanf(fid, '%d\n', 1);

s = fgets(fid);  % read a whole line

for ii = 1:nt
	for jj = 1:(NL-1)
		d(ii,jj) = fscanf(fid, '%d\t', 1);
	end
	d(ii,NL) = fscanf(fid, '%f\n', 1);
end

fclose(fid);

return;
