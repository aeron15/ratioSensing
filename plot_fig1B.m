function plot_fig1B()

%PLOT_FIG1B plots the heatmap of figrure 1 and the decision. Inspired by
%compare_thresholds_fig1B

%% LOAD settings for the comparison of S288C

gal_conc=-9:0.5:2;
glc_conc=-9:0.5:0;
gal_final = [0 2.^[gal_conc]];
glc_final = [0 2.^[glc_conc]];

gal_m = 2.^[gal_conc(1)-1:0.5:gal_conc(end)+0.5];
glc_m = 2.^[glc_conc(1)-1:0.5:glc_conc(end)+0.5];

load('../data/20140701_stitched_areas/output/M_D_stitched.mat');

%%
color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
file_names={'area','perc','mean'};

fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

% gal_m = 2.^[-10:0.5:2.5];
% glc_m = 2.^[-10:0.5:0.5];

x_tick = log2(gal_m);
y_tick = log2(glc_m);


perc_vec = [0.25];
%perc_vec = [0.2 0.35 0.5 0.65 0.8 0.95]

for i = [1 2 3]
    
    figure(i)
    h = pcolor(log2(gal_m),log2(glc_m),M{i});hold on;
    axis_handle = gca;
    set(h,'edgecolor','none');
    
    axis square
    set(axis_handle,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(axis_handle,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
    colormap(cmap)
    
    for j = 1:length(perc_vec)
        
        figure(i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        %Fit the data at the bottom by inverting the axes
        fit_cuttoff = [2^-7.5 2^-9];
        [~,~,s2] = Fit_Threshold(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        %% Parse data for plotting
        x1=log2(x);
        x1=[x1(1):0.1: x1(end)];
        y1=s(x1);
        idx_to_remove=y1<s2.b+s.b-0.1|y1>-0.5;
        y1(idx_to_remove)=[];
        x1(idx_to_remove)=[];
        
        %Plot ratio
        plot(x1,y1,'color',[0 0 0],'linewidth',2); hold on;
        %Plot threshols
        plot([s2.b s2.b],[-9 s2.b+s.b],'color',[0 0 0],'linewidth',2)
        
        axis square
        xlim([-10 2.5]);ylim([-10 0.5]);
        %set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});

    end
    
    Set_fig_RE(figure(i),17,12,12)
    title(file_names{i})
    filename=['heatmap_decision_front_WT_FIG1B_' file_names{i}];
    export_fig_specific_path(filename, '-pdf','-nocrop');
    
end

close all;
end