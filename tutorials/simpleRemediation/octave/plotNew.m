%Plot data
clear all;
close all;
nprocs = 1;
%filename = "../breakthrough.dat";
filenameb = "breakthrough_1.dat";
filenamec = "breakthrough_15.dat";
# empty separator means 'automatic'
separator = '';
skipped_rows = 0;
skipped_columns = 0;
%m = dlmread(filename, separator, skipped_rows, skipped_columns);
mb = dlmread(filenameb, separator, skipped_rows, skipped_columns);
mc = dlmread(filenamec, separator, skipped_rows, skipped_columns);
h=figure(1);
hold on;
error=0
for i=1:length(mc(:,1))
  error= error + (mb(i,2)-mc(i,2))/mc(i,2);
endfor

error/length(mc(:,1))

%loglog(m(:,1), m(:,2)/m(1,2),'-r','linewidth',3)  
loglog(mb(:,1), mb(:,2)/mb(1,2),'-b','linewidth',3)  
%loglog(mc(:,1), mc(:,2)/mc(1,2),'-r','linewidth',3)  
ylim([1e-4 1])
%xlim([50 3000])
hold off;
xlabel('Time (seconds)');
ylabel('Relative mass flow');
box on;
hold off
yt = get(gca, 'YTick');
ytkvct = 10.^linspace(1, 10*size(yt,2), 10*size(yt,2));
set(gca, 'XTick', ytkvct);
set(gca, 'XMinorTick','on', 'YMinorGrid','off')



W = 4; H = 3;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',12);

%print(h,'-deps','-color','breakthrough.eps')