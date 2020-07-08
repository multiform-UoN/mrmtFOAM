clear all;
close all;

fid = fopen("regions",'w');
l = 1;
w = 0.4;

nRegx = 19;
nRegy = 9;

kmin = 1e-15;
kmax = 1e-13;

betamin = 0.01;
betamax = 0.5;

totReg = nRegx*nRegy;

lreg = l/nRegx;
wreg = w/nRegy;

for i=1:nRegx
  
  xmin = (i-1)*lreg;
  xmax = i*lreg;
  
  for j=1:nRegy
    
    ymin = (j-1)*wreg;
    ymax = j*wreg;
    
    perm = rand*(kmax-kmin) + kmin;
    beta = rand*(betamax - betamin) + betamin;
    beta1 = rand*(betamax - betamin) + betamin;
    beta2 = rand*(betamax - betamin) + betamin;
    beta3 = rand*(betamax - betamin) + betamin;
    beta4 = rand*(betamax - betamin) + betamin;

    
    fprintf(fid,'\nboxToCell\n');
    fprintf(fid,'{\n');
    fprintf(fid,'   box (%f %f -1) (%f %f 1);\n',xmin,ymin,xmax,ymax);
    fprintf(fid,'   fieldValues\n');
    fprintf(fid,'   (\n');
    fprintf(fid,'       volScalarFieldValue beta  %f\n', beta);
    fprintf(fid,'       volScalarFieldValue beta.Sphere1  %f\n', beta1);
    fprintf(fid,'       volScalarFieldValue beta.Sphere2  %f\n', beta2);
    fprintf(fid,'       volScalarFieldValue beta.Sphere3  %f\n', beta3);
    fprintf(fid,'       volScalarFieldValue beta.Layer1  %f\n', beta4);
    fprintf(fid,'       volTensorFieldValue K (%g 0 0 0 %g 0 0 0 %g )\n',perm,perm,perm);
    fprintf(fid,'   );\n');    
    fprintf(fid,'}\n');    
    
  endfor
endfor

fflush(fid);
fclose(fid);