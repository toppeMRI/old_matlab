addpath ~/github/toppe/matlab/lib/v1

% generate files for GE scanner
seq2ge('external.seq');
system('tar czf GEscan.tgz modules.txt scanloop.txt *.mod');
%system('rm *.mod modules.txt scanloop.txt');

% simulate GE scan 
system('tar xzf GEscan.tgz');
fprintf(1, 'Simulating... ');
ny = 256;
nBlocksPerTR = 6;
d = readloop('scanloop.txt');
fprintf(1,'\n');
for ii = 1:20:ny
	fprintf(1, '\r%d of %d', ii, ny);
	scansim(1+(ii-1)*nBlocksPerTR, nBlocksPerTR+(ii-1)*nBlocksPerTR, d);
	pause(0.1);  % to allow display to update
end
fprintf(1,'\n');

