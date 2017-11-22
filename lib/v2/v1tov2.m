function v1tov2(modfile)

[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmodv1(modfile);
nomflip = round(paramsfloat(11));
mat2mod(rho,theta,gx,gy,gz,nomflip,sprintf('%sv2',modfile),desc,[],paramsint16);
