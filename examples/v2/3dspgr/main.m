writeloop;
system('./tarit');

return;

% troubleshoot THETA channel
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmod('tipdown.mod');
rho = real(rho.*exp(1i*theta));
theta = 0*theta;
mat2mod(rho,theta,gx,gy,gz,90,'tipdown,rhoonly.mod',desc);

