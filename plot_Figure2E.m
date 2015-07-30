function plot_Figure2E()

load('../data/20140617_gal80mig1_mig1mig2/output/plates_hists')


%% Mig1-Mig2-Delete

number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end

th_const = 2.5;
off_peak = 2;


gal_final = [0 2.^[-8:1:2]];
glc_final = [0 2.^[-6:1:0]];

cmap=cbrewer('seq', 'YlOrRd', 25);

%% Mig2gal80D WT


plates_hists.mig1_mig2.G12=plates_hists.mig1_mig2.G11;
plates_hists.mig1_mig2.H12=plates_hists.mig1_mig2.H11;


i=1;
glc = 2.^([-Inf   -7.5:0.5:0]);

%d = [0,0;0,12;8,0;8,12;16,0;16,12];


d=[0,0;0,0;8,8;8,8;16,16;16,16];

data = struct2cell(plates_hists);

plates={'gal80del_mig1del','mig1_mig2'};


[E_area_wt,E_prec_wt,E_mean_wt] = Plates2mat2(plates,data,plates_hists,d,map,th_const,off_peak);

%% Mig2gal80D Mt

for i=1:2
    [E_area,E_prec,E_mean] = Plates2matOther2(plates,data,plates_hists,d,map,th_const,off_peak);
    
    [E_area_wt,E_prec_wt,E_mean_wt] = Plates2mat2(plates,data,plates_hists,d,map,th_const,off_peak);

    [D_area{i},M_area{i}] = ParseHeatmapMat2(E_area{i});
    [D_mean{i},M_mean{i}] = ParseHeatmapMat2(E_mean{i});
    
    if i==1% replace area in gal80mig1
        
        M_area{1}(5,:)=M_area{1}(4,:);
        
    end
        
    figure(i)
    subplot(1,2,1)
    h = pcolor(M_area{i});
    set(h,'edgecolor','none');
    caxis([0 1])
    axis square
    set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
    set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
    colorbar
    
    subplot(1,2,2)
    h = pcolor(M_mean{i});
    %set(h,'edgecolor','none');
    axis square
    set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-8:2:2)]);
    set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-6:2:0)]);
    colorbar
    colormap(cmap)
    
    Set_fig_RE(figure(i),12,12,12);
end

end
