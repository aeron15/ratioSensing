function plot_fig1B_histograms()
%PLOT_fig1B_histograms
%cells in the population show induction (IF >0.2). 
%(Right) Comparison of cells monitored by live mi- croscopy to FCM at three sugar mixtures, denoted by numbered squares

%% Pull out specific histogram of a well in a double gradient
% from replicate 2 (-6.5,-8) (-6.5,-2)

p = 1; i=7;j = 4;

plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
load('../data/s288c/20140217_4283/output/plates_hists_EMD')

data = struct2cell(plates_hists);

plate_name = 'Plate_A';well_name = 'G04';
[x,y] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_A';well_name = 'G10';
[x2,y2] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_C';well_name = 'D10';
[x3,y3] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_D';well_name = 'D02';
[x4,y4] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

max_x = max([max(x),max(x2),max(x3),max(x4)]);
min_x = min([min(x),min(x2),min(x3),min(x4)]);
xx = linspace(min_x,max_x,150);

yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
yy(isnan(yy))=0;

yy2 = interp1(x2,y2,xx,'linear');yy2(yy2<0)=0;
yy2(isnan(yy2))=0;

yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
yy(isnan(yy))=0;

yy2 = interp1(x2,y2,xx,'linear');yy2(yy2<0)=0;
yy2(isnan(yy2))=0;

yy3 = interp1(x3,y3,xx,'linear');yy3(yy3<0)=0;
yy3(isnan(yy3))=0;

yy4 = interp1(x4,y4,xx,'linear');yy4(yy4<0)=0;
yy4(isnan(yy4))=0;

%%
figure(10)
plot(xx,yy/sum(yy),'k','linewidth',6);box off;
%ylim([0 0.08])
Set_fig_RE(figure(10),25,18,18)
filename='Figure1B_hist_1';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');


figure(11)
plot(xx,yy2/sum(yy2),'k','linewidth',6);box off;
Set_fig_RE(figure(11),25,18,18)
%ylim([0 0.08])
filename='Figure1B_hist_2';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');


figure(12)
plot(xx,yy3/sum(yy3),'k','linewidth',6);box off;
Set_fig_RE(figure(12),25,18,18)
%ylim([0 0.08])
filename='Figure1B_hist_3';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

% figure(13)
% plot(xx,yy4,'k','linewidth',6);box off;
% Set_fig_RE(figure(13),25,18,18)

end