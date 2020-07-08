%Plot data
clear all;
close all;
nprocs = 1;
filename = "breakthrough_spheresLayers.dat";
filenameb = "breakthrough_spheres.dat";
filenamec = "breakthrough_no.dat";
# empty separator means 'automatic'
separator = '';
skipped_rows = 0;
skipped_columns = 0;
m = dlmread(filename, separator, skipped_rows, skipped_columns);
mb = dlmread(filenameb, separator, skipped_rows, skipped_columns);
mc = dlmread(filenamec, separator, skipped_rows, skipped_columns);
ttot = m(:,1);
etatot_ = m(:, 2);


ttotb = mb(:,1);
etatotb_ = mb(:, 2);

for i=1:floor(length(etatot_))/nprocs

  t(i) = ttot(i*nprocs);
  eta_(i) = etatot_(i*nprocs);
  

endfor

for i=1:floor(length(etatotb_))/nprocs

  tb(i) = ttotb(i*nprocs);
  etab_(i) = etatotb_(i*nprocs);
  

endfor
h=figure(1);
hold on;

loglog(t, eta_/eta_(1),'-r','linewidth',3)  
loglog(tb, etab_/etab_(1),'-b','linewidth',3)  
loglog(mc(1:500,1), mc(1:500,2)/mc(1,2),'--k','linewidth',3)  
%semilogy(tb, etab_/etab_(1),'-b')  

%semilogy(m2(:,1)/50, m2(:,2),'*r') 
%semilogy(m3(:,1), m3(:,2),'*b')  
xlim ([0, 550])
ylim ([0.001,1])
xlabel('Time (seconds)');
ylabel('Relative mass flow');
box on;
hold off
L = legend("Spheres + Layer ",...
          "Spheres ",...
          "No MRMT")
legend boxoff;
FL1= findall(L,'-property','FontName');
FL2= findall(L,'-property','FontSize');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(FL2,'FontSize',12);
W = 4; H = 3;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',12);

print(h,'-deps','-color','breakthrough.eps')