% addpath /net/brooks/export/home/jfnielse/projects/rfovfmri/b0map/recon/\" >>$bindir/$scriptname.m"); system($systemstr);
addpath /net/brooks/export/home/jfnielse/github/toppe/matlab/lib/GE/  % load raw data from Pfile
addpath /net/brooks/export/home/jfnielse/github/toppe/matlab/lib/v2/  % readmod
curdir = pwd; cd /net/brooks/export/home/jfnielse/projects/autorfdesign/SPINS2/mtools/irt; setup; cd(curdir);

