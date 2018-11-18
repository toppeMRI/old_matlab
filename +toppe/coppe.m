% Packages toppe files into toppe-scanfiles.tgz, then copies it to the scanner
% Assumes you have SSH keys set up to log into romero
disp('Making archive...');
system('tar czf toppe-scanfiles.tgz modules.txt scanloop.txt *.mod');
disp('Copying to romero...');
system('scp toppe-scanfiles.tgz fmrilab@romero:~/amos/');
disp('Copying to scanner...');
system('ssh -q fmrilab@romero "~/amos/pushtoppefiles"');
disp('Done!');