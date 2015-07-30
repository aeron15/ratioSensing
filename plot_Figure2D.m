function plot_Figure2D()
%PLOT_FIGURE2D
%n a mig1? background, cells continue to respond to the galactose:glucose ratio. 
%Experiment performed in duplicate. Solid line represents the decision front of the mig1?; 
%dashed line represents the decision front of the wild-type strain 

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

load('../data/20140508_mig1_del/output/plates_hists')

%% mig1  WT
glc = [2.^[-8:1:1]];
gal = [2.^[-10:1:3]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'mig1_del'};

[E_area{1},E_prec{1},E_mean{1}] = Plates2matOther(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

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
filename=['Figure_2C_mig1_delete'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');
%%
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

% Compare mig1D to WT
glc = [2.^[-6:1:0]];
gal = [2.^[-8:1:2]];

fit_cuttoff = [2^-1 2^-9];
mid_value = 2^-4;
cutoff=0.2;
color_vec = [0 0 0;0 0.5 0]

clear a b a_d a_u b_d b_u
%% WT
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
Set_fig_RE(figure(1),17,12,12)


%%
close all;

end