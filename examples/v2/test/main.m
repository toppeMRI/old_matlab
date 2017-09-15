
addpath ~/github/toppe/matlab/lib/v1
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod('tipdown_real.mod');
rmpath ~/github/toppe/matlab/lib/v1
addpath ~/github/toppe/matlab/lib/v2
mat2mod(rho,0*rho,gx,gy,gz,25,'tipdown_real_v2.mod',desc);
mat2mod(abs(rho),angle(rho),gx,gy,gz,25,'tipdown_abs_v2.mod',desc);

writeloop;
system('./tarit');
