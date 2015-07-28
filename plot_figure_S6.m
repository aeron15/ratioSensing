function plot_figure_S6()

%% LOAD settings for the comparison of S288C

gal_conc=-9:0.5:2;
glc_conc=-9:0.5:0;
gal_final = [0 2.^[gal_conc]];
glc_final = [0 2.^[glc_conc]];

gal_m = 2.^[gal_conc(1)-1:0.5:gal_conc(end)+0.5];
glc_m = 2.^[glc_conc(1)-1:0.5:glc_conc(end)+0.5];

load('../data/20140701_stitched_areas/output/M_for_replicate_original_supp_material.mat')
load('../data/20140701_stitched_areas/output/M_D_stitched.mat')

%% Plot figure 

color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
file_names={'area','perc','mean'};

fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

x_tick = log2(gal_m);
y_tick = log2(glc_m);

perc_vec = [0.2 0.35 0.5 0.65 0.8 0.95]

for i = [1 2 3]
    
    figure(i)
    h = pcolor(log2(gal_m),log2(glc_m),M{i});hold on;
    axis_handle = gca;
    set(h,'edgecolor','none');
    
    axis square
    set(axis_handle,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(axis_handle,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
    h = colorbar('location','northoutside');
    colormap(cmap)
    Set_fig_RE(figure(i),17,12,12)
    
    title(file_names{i})
    filename=['Heatmap_S6_' file_names{i}];
    box off;
    export_fig_specific_path(filename, '-pdf','-nocrop');
    
    for j = 1:length(perc_vec)
        
        figure(10*i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        plot(log2(x),log2(y),'o','markerfacecolor',color_vec(j,:),'color',color_vec(j,:),'markersize',4);hold on
        axis square
        xlim([-10 2.5]);ylim([-10 0.5]);
        
        set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
        Set_fig_RE(figure(10*i),17,12,12)
        
    end
    
    title(file_names{i})
    filename=['contour_WT_S6_' file_names{i}];
    box off;
    export_fig_specific_path(filename, '-pdf','-nocrop');
    close all;
    
end