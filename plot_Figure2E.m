function plot_Figure2E()

load('../data/20140617_gal80mig1_mig1mig2/output/plates_hists')

%% gal80mig1 and mig1mig2

map =create_map;

th_const = 2.5;
off_peak = 2;

gal_final = [0 2.^[-8:1:2]];
glc_final = [0 2.^[-6:1:0]];

cmap=cbrewer('seq', 'YlOrRd', 25);

plates_hists.mig1_mig2.G12=plates_hists.mig1_mig2.G11;
plates_hists.mig1_mig2.H12=plates_hists.mig1_mig2.H11;

glc = 2.^([-Inf   -7.5:0.5:0]);

d=[0,0;0,0;8,8;8,8;16,16;16,16];

data = struct2cell(plates_hists);

plates={'gal80del_mig1del','mig1_mig2'};

[E_area_wt,E_prec_wt,E_mean_wt] = Plates2mat2(plates,data,plates_hists,d,map,th_const,off_peak);

%% Mig1gal80D  and mig1mig2D

for i=1:2
    [E_area,E_prec,E_mean] = Plates2matOther2(plates,data,plates_hists,d,map,th_const,off_peak);
    
    [E_area_wt,E_prec_wt,E_mean_wt] = Plates2mat2(plates,data,plates_hists,d,map,th_const,off_peak);
    
    [D_area{i},M_area{i}] = ParseHeatmapMat2(E_area{i});
    [D_mean{i},M_mean{i}] = ParseHeatmapMat2(E_mean{i});
    
    if i==1% replace area in gal80mig1
        M_area{1}(5,:)=M_area{1}(4,:);
    end
    
    figure(i)
    h = pcolor(M_area{i});
    set(h,'edgecolor','none');
    caxis([0 1])
    axis square
    set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
    set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
    colorbar
    colormap(cmap)
    
    switch i
        
        case 1
            title('\it{mig1\Delta gal80\Delta}')
            filename = 'Figure_2E_mig1_gal80';
            
        case 2
            title('\it{mig1\Delta mig2\Delta}')
            filename = 'Figure_S5C_mig1_mig2';
    end
    
    Set_fig_RE(figure(i),12,12,12);
    export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');
    
end

close all;

end
