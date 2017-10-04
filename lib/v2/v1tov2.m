function v1tov2(modfile)

[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmodv1(modfile);
mat2mod(rho,theta,gx,gy,gz,90,sprintf('%sv2',modfile),desc,[],paramsint16);
