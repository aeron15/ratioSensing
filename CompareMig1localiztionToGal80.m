function CompareMig1localiztionToGal80

%% plot gal80 compare to mig1 localization

number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end

th_const = 2.5;
off_peak = 2;


gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-9:0.5:0]];

cmap=cbrewer('seq', 'YlOrRd', 25);

%%

%% Gal80  delete
% close all

load('../data/20140407_gal80_del/output/plates_hists')

glc = [2.^[-8:1:1]];
gal = [2.^[-10:1:3]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'Gal80_del'};



% GAl2 rep1
i = 1;

[E_area{i},E_prec{i},E_mean{i}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);



[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(5)
h = pcolor(log2(gal),log2(glc),M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+(-11)+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 3 5:2:22]+(-9)+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
title('\it{gal80\Delta} Area')
colormap(cmap)
caxis([0, 1])
title('\it{gal80\Delta}')

Set_fig_RE(figure(5),12,12,12)

%% Compute the average across the rows

average_mean=nanmean(E_mean{i}(2:end,:),2);

average_plotted=[average_mean average_mean];

average_plotted=[average_plotted; average_plotted(end,:)];

figure(7)
h = pcolor(average_plotted);

set(h,'edgecolor','none');
%axis square
colormap(cmap)
%caxis([0, 1])






%%

figure(6)
h = pcolor(log2(gal),log2(glc),M_mean{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+(-11)+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 3 5:2:22]+(-9)+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
title('\it{gal80\Delta} Area')
colormap(cmap)
% caxis([0, 1])
title('\it{gal80\Delta}')

Set_fig_RE(figure(6),12,12,12)

%% Compute the average across the rows

average_area=nanmean(E_area{i}(2:end,:),2);

average_plotted=[average_area average_area];

average_plotted=[average_plotted; average_plotted(end,:)];

hfig=figure('Position',[440   378   103   420])
h = pcolor(average_plotted);

set(h,'edgecolor','none');
%axis square
colormap(cmap)
caxis([0, 1])
Set_fig_RE(hfig,18,18,22);

filename='Average_induction_gal80del_for_glucose_concentrations';
export_fig(filename, '-pdf','-transparent','-nocrop');

%% Compare to Mgf1p localization. Done by automated segmentation and also by hand.

 ratio_M =  [ 0.7489    0.8010    1.0000    0.8885    0.7418    0.9128
              0.7450    0.6293    0.7638    0.8655    0.8481    0.7859
              0.1976    0.4528    0.2202    0.1088    0.0285    0.0920
              0.2550    0.1662    0.2701    0.2662         0    0.1539];

          
          
         

   ratio_M_count =  [1.0000    0.9231    0.7059    0.7857    0.6842    0.8919
    0.7097    0.6452    0.7143    0.6667    0.8000    0.8529
    0.0952    0.1667    0.0769    0.3000    0.0909    0.0690
         0         0         0         0         0    0.1000];

         
figure(1);hold on;
glc_80 = [ -6:0 ];       
plot(glc_80,mean(mat2gray(D_mean{1}),2),'-ro','markerface','r')
      
% plot(glc_80,mean(mat2gray(D_area{1}),2),'ro','markerface','g')


glc_mig1 = [-8 -7 -5 -3]
% plot(glc_mig1,mean(1-ratio_M,2),'ok','markerface','k')
plot(glc_mig1,mean(1-ratio_M_count,2),'-ok','markerface','k')
glc = [-7:0];
set(gca,'xtick',[-8:0],'xticklabel',['No|',sprintf('%0.1g|', (2.^(glc)))]);
xlabel('Glucose [%]');
ylabel('Response');
%legend({'GAL 1 Mean Expression (\itgal80\Delta\rm)','Fraction of localized cells'});
Set_fig_RE(figure(1),12,12,18)

filename=[' Figure_Mig1p_gal80_delete_response_2'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%print(figure(1),'-dtiff','-r300')

%%

% reposne on the x-axis
figure(2);hold on;
glc_80 = [ -6:0 ];    
response_80 = mean(mat2gray(D_mean{1}),2);

plot(response_80(:),glc_80(:),'-ro','markerface','r')

glc_mig1 = [-8 -7 -5 -3]
% plot(glc_mig1,mean(1-ratio_M,2),'ok','markerface','k')
plot(mean(1-ratio_M_count,2),glc_mig1,'-ok','markerface','k');

set(gca,'ytick',[-8:0],'yticklabel',['No|',sprintf('%0.1g|', (2.^(glc)))]);


glc = [-7:0];
ylabel('Glucose [%]');
xlabel('Response');
legend({'GAL 1 Mean Expression (\itgal80\Delta\rm)','Fraction of localized cells'});
Set_fig_RE(figure(2),12,12,18)

filename=['Figure_Mig1p_gal80_delete_response'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%print(figure(2),'-dtiff','-r300')
