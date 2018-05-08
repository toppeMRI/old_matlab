addpath ../../v2

% generate files for GE scanner
seq2ge('external.seq');
system('tar czf GEscan.tgz modules.txt scanloop.txt *.mod');
%system('rm *.mod modules.txt scanloop.txt');

% simulate GE scan 
nCoresPerTR = 6;
tpause = 0.05;    % sec (slow down display loop)
nBlocksPerTR = 6;
playseq('scanloop.txt',nBlocksPerTR,0,tpause);

return;

