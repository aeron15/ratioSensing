function figure3_YS

%% Plots for figure3. Gal80p and Mig1p deleteions

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
% cmap=jet(100);

% path_figures='/Users/RenanEscalante/Dropbox/Galactose_Pathway/gal_paper/figures/Figure_2/';

%% mig1  delete other
close all
%Yoni
%load('C:\Users\Yoni\Dropbox\gal_paper\Data\20140508_mig1_del\output\plates_hists')
load('../data/20140508_mig1_del/output/plates_hists')

%Renan
% load('/Users/RenanEscalante/Dropbox/Galactose_Pathway/gal_paper/Data/20140508_mig1_del/output/plates_hists')

glc = [2.^[-8:1:1]];
gal = [2.^[-10:1:3]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'mig1_del'};


[E_area{1},E_prec{1},E_mean{1}] = Plates2matOther(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

%
i = 1;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(log2(gal),log2(glc),M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+(-11)+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 3 5:2:22]+(-9)+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
title('\it{mig1\Delta} Area')

% colorbar 
colormap(cmap)
caxis([0, 1])

Set_fig_RE(figure(1),12,12,12)
filename=['Figure_2C_mig1_delete_jet'];
% export_fig(filename, '-pdf','-transparent','-nocrop');

i = 1;
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(2)
h = pcolor(M_mean{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-6:2:0)]);
title('\it{mig1\Delta} Mean')

colorbar 
colormap(cmap)

Set_fig_RE(figure(2),12,12,12)

% WT


[E_area{2},E_prec{2},E_mean{2}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

%
i = 2;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(3)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-6:2:0)]);

colorbar 
colormap(cmap)

Set_fig_RE(figure(3),12,12,12)
i = 2;
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(4)
h = pcolor(M_mean{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-6:2:0)]);

colorbar 
colormap(cmap)

Set_fig_RE(figure(4),12,12,12)

% Compare mig1D to WT
glc = [2.^[-6:1:0]];
gal = [2.^[-8:1:2]];

fit_cuttoff = [2^-1 2^-9];
mid_value = 2^-4;
cutoff=0.2;
color_vec = [0 0 0;0 0.5 0]

clear a b a_d a_u b_d b_u
% WT
figure(1);hold on
i=1
for j = [1 2]
   [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{j},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
%     plot(log2(x),log2(y),'o','markerfacecolor',color_vec(i,:),'color',color_vec(i,:),'markersize',2);hold on
x = gal; y = s(log2(x)); ind  = find(y>=(-6));
    plot(log2(x(ind)),y(ind),'color',color_vec(i,:),'linewidth',4); hold on;
    xlim([-10 3]);ylim([-8 1]);
    i = i+1;
end
% figure(4);axis square
Set_fig_RE(figure(4),17,12,12)

% filename=['Figure_MIG1'];
% export_fig(filename, '-pdf','-transparent','-nocrop');



%% Gal80  delete
% close all
%load('C:\Users\Yoni\Dropbox\gal_paper\Data\20140407_gal80_del\output\plates_hists')
load('../data/20140407_gal80_del/output/plates_hists')

%load('/Users/RenanEscalante/Dropbox/Galactose_Pathway/gal_paper/Data/20140407_gal80_del/output/plates_hists')

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

% filename=['Figure_2C_gal80_delete_jet'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
% gal80D WT
% load('C:\Users\Yoni\Dropbox\gal_paper\Data\20140407_gal80_del\output\plates_hists')
%Renan
%load('../../Data/20140407_gal80_del/output/plates_hists')
% load('../../../Data/20140407_gal80_del/output/plates_hists')
d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'Gal80_del'};

i = 2;
[E_area{i},E_prec{i},E_mean{i}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

%[E_area{i},E_prec{i},E_mean{i}] = Plates2matBFP(plates,data,plates_hists,d,map,th_const,off_peak);%modified by Rena to for a not BFP version

[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(6)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
title('WT grown with \it{gal80\Delta}')

colorbar 
colormap(cmap)
caxis([0, 1])

Set_fig_RE(figure(6),12,12,12)

% filename=[path_figures 'Figure_2C_gal80_delete_WT'];
% export_fig(filename, '-pdf','-transparent','-nocrop');


% Compare gal80D to WT
glc = [2.^[-6:1:0]];
gal = [2.^[-8:1:2]];

fit_cuttoff = [2^-1 2^-9];
mid_value = 2^-4;
cutoff=0.2;
color_vec = [0 0 0;0 0.5 0]

clear a b a_d a_u b_d b_u
% WT
figure(5);hold on
i=1
for j = [1 2]
   [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{j},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
%     plot(log2(x),log2(y),'o','markerfacecolor',color_vec(i,:),'color',color_vec(i,:),'markersize',2);hold on
    if j==1;
        s = fit(log2(x),log2(y),'a*x+b');
    end
    
x = gal; y = s(log2(x)); ind  = find(y>=(-6));
    plot(log2(x(ind)),y(ind),'color',color_vec(i,:),'linewidth',4); hold on;
    xlim([-10 3]);ylim([-8 1]);
    i = i+1;
end
% figure(4);axis square
% Set_fig_RE(figure(4),17,12,12)

%% gal80 counts

home = pwd;
cd ('C:\Users\Yoni\Dropbox\Mig1\matlab');

figure(7)
i = 1
[E_area{i},E_prec{i},E_mean{i},E_mean_area{i},E_counts{i}] = Plates2matMchSingle(plates{1},data(1),plates_hists,d,map,th_const,off_peak);

[D_counts{i},M_counts{i}] = ParseHeatmapMat(E_counts{i});

h = pcolor(M_counts{i});
set(h,'edgecolor','none');
axis square

set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No|',sprintf('%0.1g|', (2.^(-8:2:2)))]);
set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No|',sprintf('%0.1g|',2.^(-6:2:0))]);




colorbar 
colormap(jet)
% caxis([0, 1])
figure(7);axis square
Set_fig_YS(figure(7),17,12,12)

%% GAl2 rep1
% 
% load('C:\Users\Yoni\Dropbox\gal_paper\Data\20140407_gal80_del2\output\plates_hists')
% % load('/Users/RenanEscalante/Dropbox/Galactose_Pathway/gal_paper/Data/20140407_gal80_del2/output/plates_hists')
% d = [0,0;0,12;8,0;8,12;16,0;16,12];
% 
% data = struct2cell(plates_hists);
% plates = {'Gal80_del_2'};
% 
% i=2
% [E_area{i},E_prec{i},E_mean{i}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);
% 
% % remove low and high gloucse rows:
% [D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
% 
% figure(i)
% h = pcolor(M_area{i});
% set(h,'edgecolor','none');
% axis square
% set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-8:2:2)]);
% set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-6:2:0)]);
% 
% colorbar 
% colormap(cmap)
% 
% Set_fig_RE(figure(i),12,12,12)



%%
% 
% figure(7)
% plot(D_area_wt{1}','k');hold on;plot(D_area{1}','r');
% figure(8)
% plot(D_mean_wt{1}','k');hold on;plot(D_mean{1}','r');
% 
% %% Gal80  delete mean
% 
% i=1
% [D_area{i},M_area{i}] = ParseHeatmapMat(E_mean{i});
% 
% figure(3)
% h = pcolor(M_area{i});
% set(h,'edgecolor','none');
% axis square
% set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
% set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
% 
% colorbar 
% colormap(cmap)
% 
% Set_fig_RE(figure(1),12,12,12)
% 
% 
% 
% i=2
% 
% [D_area{i},M_area{i}] = ParseHeatmapMat(E_mean{i});
% 
% figure(4)
% h = pcolor(M_area{i});
% set(h,'edgecolor','none');
% axis square
% set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-8:2:2)]);
% set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-6:2:0)]);
% 
% colorbar 
% colormap(cmap)
% 
% Set_fig_RE(figure(i),12,12,12)
% 
% %%
% figure(5)
% plot(D_area{1}')
% 
% 
% %% Analyze the GEV strain
% color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
% 
% close all
% %load('C:\Users\ys151\Dropbox\gal_paper\Data\20140520_gal2del_GEV\output\plates_other')
% load('/Users/RenanEscalante/Dropbox/Galactose_Pathway/gal_paper/Data/20140520_gal2del_GEV/output/plates_other')
% 
% ctrs = logspace(-4,6,200);
% estordial = ([128 64 32 16 8 4 2 0]);
% d = [0,0;0,7;4,0;4,7];
% 
% 
% data = struct2cell(plates_other);
% xx = struct2cell(data{1});
% 
% for i =1:length(xx)
% 
%     y_serial(i,:) = hist(xx{i}.yfp,ctrs);
%     E_y_serial(i) = mean(xx{i}.yfp);
%     std_y_serial(i) = std(xx{i}.yfp);
% 
% end
% 
% figure(1)
% plot(log10(ctrs),y_serial);xlim([0 5])
% errorbar(estordial,E_y_serial,std_y_serial);set(gca,'yscale','log');
% ylabel('log10(Gal1 YFP)');xlabel('E');
% Set_fig_RE(figure(1),12,12,12);
% 
% plates = {'small_DG'};
% 
% data = struct2cell(plates_other);
% [M_hist,M_mean{3},M_std] = Plates2matGEV(plates,data(2),plates_other,d,map,th_const,off_peak);
% 
% x = flipud(padarray(log10(M_mean{3}),[1 1],nan));
% figure(2);pcolor(x(2:end,2:end));
% axis square;colormap(cmap);colorbar;
% Set_fig_RE(figure(2),12,12,12);
% 
% glc_gev = fliplr([-6 -4 -2 0]);
% glc_80 = [-7:-1];
% 
% 
% for b = 1:4
%     figure(3)
%     subplot(2,2,b)
%         plot(glc_80,D_mean{1},'o-k');hold on;
%     plot(glc_gev,log10(M_mean{3}([1:4]+d(b,1),[1:4]+d(b,2))),'o-','color',color_vec(b,:));hold on;
% 
%     ylabel('log10(Gal1 YFP)');xlabel('Glc');
% 
%     ylim([1.5 4.5]);xlim([-7 0]);
% end
% Set_fig_RE(figure(3),12,12,12);
% 
% color_vec = [1 0 0; 0 1 0; 0 0 1;1 1/2 0]
% for b = 1:4
%     figure(4)
%         plot(glc_80,D_mean{1},'o-k');hold on;
%     plot(glc_gev,log10(M_mean{3}([1:4]+d(b,1),[1:4]+d(b,2))),'o-','color',color_vec(b,:));hold on;
%     ylabel('log10(Gal1 YFP)');xlabel('Glc');
% 
%     ylim([1.5 4.5]);xlim([-7 0]);
% end
% Set_fig_RE(figure(4),12,12,12);