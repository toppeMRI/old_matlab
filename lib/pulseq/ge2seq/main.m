addpath ~jfnielse/github/toppe/matlab/lib     % readmod.m
addpath ~jfnielse/projects/pulseq/lib/matlab/         % Pulseq Matlab toolbox

ge2seq('modules.txt', 'scanloop.txt', 'timing.txt');
