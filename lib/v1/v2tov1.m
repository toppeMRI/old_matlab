function v2tov1(modfile)

[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmodv2(modfile);
nomflip = round(paramsfloat(11));
mat2mod(rho,theta,gx,gy,gz,nomflip,sprintf('%sv1',modfile),desc,[],paramsint16);
