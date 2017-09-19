% convert .mod files from v1 to v2
rmpath ~jfnielse/github/toppe/matlab/lib/v2/
addpath ~jfnielse/github/toppe/matlab/lib/v1/
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod('tipdown,v1.mod');
rmpath ~jfnielse/github/toppe/matlab/lib/v1/
addpath ~jfnielse/github/toppe/matlab/lib/v2/
mat2mod(rho,theta,gx,gy,gz,90,'tipdown.mod',desc);

rmpath ~jfnielse/github/toppe/matlab/lib/v2/
addpath ~jfnielse/github/toppe/matlab/lib/v1/
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod('rf180,v1.mod');
rmpath ~jfnielse/github/toppe/matlab/lib/v1/
addpath ~jfnielse/github/toppe/matlab/lib/v2/
mat2mod(rho,theta,gx,gy,gz,180,'rf180.mod',desc);

% Create scan.tgz. Untar to /usr/g/bin/ on scanner and scan with toppev2
writeloop;
system('./tarit');

% display sequence
playseq('scanloop.txt',3,0,1);
