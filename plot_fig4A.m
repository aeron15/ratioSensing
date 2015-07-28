function plot_fig4A()
%PLOT_FIG4A. A wild-type strain was competed against a gal4? in two concentrations of carbon:
%Inspired by compare_thresholds_fig_fitness_schematic

%% LOAD settings for the comparison of S288C
gal_conc=-9:0.5:2;
glc_conc=-9:0.5:0;
gal_final = [0 2.^[gal_conc]];
glc_final = [0 2.^[glc_conc]];

gal_m = 2.^[gal_conc(1)-1:0.5:gal_conc(end)+0.5];
glc_m = 2.^[glc_conc(1)-1:0.5:glc_conc(end)+0.5];

%load('../data/20140701_stitched_areas/output/M_for_replicate_original_supp_material.mat')
load('../data/20140701_stitched_areas/output/M_D_stitched.mat')

%%
color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
file_names={'area','perc','mean'};

fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;

green=[13 144 69];
gray=[51 56 53];

green=green./255;
gray=gray./255;


cmap(2,:)=green;
cmap(1,:)=gray;

x_tick = log2(gal_m);
y_tick = log2(glc_m);

perc_vec = [0.22];

for i = [1]%1 is the area metric [1 2 3] to plot also mean and percentage
    
    figure(i)
    
    area=M{i};
    area(area>-0.1)=nan;
    area(9,1)=1;
    area(9,23)=0;
    h = pcolor(log2(gal_m),log2(glc_m),area);hold on;
    axis_handle = gca;
    set(h,'edgecolor','none');
    
    axis square
    set(axis_handle,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(axis_handle,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
    colormap(cmap)
    
    for j = 1:length(perc_vec)
        
        figure(i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        x1=log2(x);
        x1=[x1(1):0.1: x1(end)];
        y1=s(x1);
        idx_to_remove=y1<-7.8|y1>-0.5;
        y1(idx_to_remove)=[];
        x1(idx_to_remove)=[];
        
        plot(x1,y1,'color',[0 0 0],'linewidth',2); hold on;
        plot([-8 -8],[-10 -7.26],'color',[0 0 0],'linewidth',2)
        plot([-11 2],[-7.26 -7.26],'color',[0 0 0],'linewidth',2)
        
    end
    
    Set_fig_RE(figure(i),18,18,18)

    xlim([-10 1.5]);
    ylim([-9 -2.5]);
    title(file_names{i})
    set(gcf, 'color', 'none');
    set(gca, 'color', 'none');
    filename=['contour_heatmap_WT_fitness' file_names{i}];
    export_fig_specific_path(filename, '-pdf','-nocrop');
    
end

close all;
end