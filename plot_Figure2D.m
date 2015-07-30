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

%% mig1  WT
load('../data/20140508_mig1_del/output/plates_hists')
gal = [0 2.^[-9:0.5:2]];
glc = [0 2.^[-10.5:0.5:0]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'mig1_del'};


[E_area{1},E_prec{1},E_mean{1}] = Plates2matOther(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:
i = 1;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_area{i},M_area{i}] = ParseHeatmapMat(E_prec{i});

figure('Name','mig1 area')
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', {-8:2:2}]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', {-6:2:0}]);

colorbar 
colormap(cmap)

Set_fig_RE(figure(1),12,12,12)
%%
filename='Figure_2D_mig1delete';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

close all;

end