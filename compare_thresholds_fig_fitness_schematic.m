function compare_thresholds_fig_fitness_schematic(gal_final,glc_final,gal_m,glc_m,M,D)

%%

color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
file_names={'area','perc','mean'};

fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
%cmap=cbrewer('seq', 'YlOrRd', 25);
green=[13 144 69];
gray=[51 56 53];

green=green./255;
gray=gray./255;

% cmap(2,:)=rgb('SeaGreen');
% cmap(1,:)=rgb('Gray');

cmap(2,:)=green;
cmap(1,:)=gray;

% gal_m = 2.^[-10:0.5:2.5];
% glc_m = 2.^[-10:0.5:0.5];

x_tick = log2(gal_m);
y_tick = log2(glc_m);


perc_vec = [0.22];
%perc_vec = [0.2 0.35 0.5 0.65 0.8 0.95]

for i = [1 2 3]
    
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
    %h = colorbar('location','northoutside');
    colormap(cmap)
    
    for j = 1:length(perc_vec)
        
        figure(i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        %plot(log2(x),log2(y),'o','markerfacecolor',[0 0 0],'color',[0 0 0],'markersize',4);hold on
        x1=log2(x);
        x1=[x1(1):0.1: x1(end)];
        y1=s(x1);
        idx_to_remove=y1<-7.8|y1>-0.5;
        y1(idx_to_remove)=[];
        x1(idx_to_remove)=[];
        
        plot(x1,y1,'color',[0 0 0],'linewidth',2); hold on;
                
        plot([-8 -8],[-10 -7.26],'color',[0 0 0],'linewidth',2)
        
        plot([-11 2],[-7.26 -7.26],'color',[0 0 0],'linewidth',2)

        %set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
        
    end
    
    Set_fig_RE(figure(i),18,18,18)

    xlim([-10 1.5]);
    ylim([-9 -2.5]);
    %axis square
    %hline(-7.26)    
    title(file_names{i})
    set(gcf, 'color', 'none');
    set(gca, 'color', 'none');
    filename=['contour_heatmap_WT_fitness' file_names{i}];
    %export_fig(filename, '-pdf','-nocrop');
    close all;
    
end

end