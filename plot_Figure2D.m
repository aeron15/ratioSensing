function plot_Figure2D()
%PLOT_FIGURE2D
%n a mig1? background, cells continue to respond to the galactose:glucose ratio.
%Experiment performed in duplicate. Solid line represents the decision front of the mig1?;
%dashed line represents the decision front of the wild-type strain
% mCherry is the wild type strain
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

load('../data/20140508_mig1_del/output/plates_hists')

%% mig1 delete and  WT
glc = [2.^[-8:1:1]];
gal = [2.^[-10:1:3]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'mig1_del'};

[E_area{1},~,~] = Plates2matOther(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high glucose rows:
i = 1;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(log2(gal),log2(glc),M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+(-11)+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
set(gca,'ytick',[1 3 5:2:22]+(-9)+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
title('\it{mig1\Delta} Area')

colormap(cmap)
caxis([0, 1])

%% Compare mig1D to WT

% WT
[E_area{2},~,~] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);

i = 2;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

glc = [2.^[-6:1:0]];
gal = [2.^[-8:1:2]];

fit_cuttoff = [2^-1 2^-9];
mid_value = 2^-4;
cutoff=0.2;
color_vec = [0 0 0;0 0.5 0]

%% Add the decision fronts of WT in green and Mig1 delete in black
% 1 is the mig1 delete and 2 is the WT strain
figure(1);hold on
i=1
for j = [1 2]
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{j},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
    x = gal; y = s(log2(x)); ind  = find(y>=(-6));
    plot(log2(x(ind)),y(ind),'color',color_vec(i,:),'linewidth',4); hold on;
    xlim([-10 3]);ylim([-8 1]);
    i = i+1;
end

%%
Set_fig_RE(figure(1),17,12,12)
filename=['Figure_2C_mig1_delete'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');
close all;

end