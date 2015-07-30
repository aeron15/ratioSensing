function plot_Figure2C
%PLOT_FIGURE2C
%In a gal80? background, the ratio response is converted to a threshold response 
%(i.e., in the absence of Gal80p the response is galactose independent). 
%Experiment performed in duplicate. Data for no glucose conditions is not shown for clarity (Methods)
% based on figure3_YS

%% Create map of 96 well plate
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
 
%% Gal80  delete

load('../data/20140407_gal80_del/output/plates_hists')

%%

glc = [2.^[-8:1:1]];
gal = [2.^[-10:1:3]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'Gal80_del'};

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
filename='Figure_2C_gal80delete';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');


%% Plot wild type strain
% d = [0,0;0,12;8,0;8,12;16,0;16,12];
% 
% data = struct2cell(plates_hists);
% plates = {'Gal80_del'};
% 
% i = 2;
% [E_area{i},E_prec{i},E_mean{i}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);
% 
% 
% [D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
% [D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});
% 
% figure(6)
% h = pcolor(M_area{i});
% set(h,'edgecolor','none');
% axis square
% set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
% set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
% title('WT grown with \it{gal80\Delta}')
% 
% colorbar 
% colormap(cmap)
% caxis([0, 1])
% 
% Set_fig_RE(figure(6),12,12,12)
% 
% filename=['Figure_2C_gal80_delete_WT'];
%export_fig(filename, '-pdf','-transparent','-nocrop');


close all;

