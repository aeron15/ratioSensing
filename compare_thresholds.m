function compare_thresholds(gal_final,glc_final,gal_m,glc_m,M,D)

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


perc_vec = [0.2 0.35 0.5 0.65 0.8 0.95]

for i = [1 2 3]
    
    %    M{i}(7,3) = mean([M{i}(6,3),M{i}(8,3)]);
    %    M_area{i}(1,1) = 1;
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
    
%     if i==3
%         
%         % To correct a small imperfection at the bottom of the mean plot
%         perc_vec = [0.2 0.35 0.5 0.65 0.8 0.88]
%         
%     end
    
    for j = 1:length(perc_vec)
        
        figure(10*i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,perc_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        
        plot(log2(x),log2(y),'o','markerfacecolor',color_vec(j,:),'color',color_vec(j,:),'markersize',4);hold on
        %plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
        axis square
        xlim([-10 2.5]);ylim([-10 0.5]);
        
        set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
        Set_fig_RE(figure(i),17,12,12)
        Set_fig_RE(figure(10*i),17,12,12)
        
    end
    %print('-dtiff','-r1000',[figure_folder,'measures',num2str(i)]);
    title(file_names{i})
    filename=['contour_WT_' file_names{i}];
    box off;
    export_fig_specific_path(filename, '-pdf','-nocrop');
    close all;
    
end

end