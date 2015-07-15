function figureSI4()

%% Figure SI4 heterozygote delete

number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end


th_const = 2.5;
off_peak = 2;

cmap=cbrewer('seq', 'YlOrRd', 25);
figure_folder = '.';
fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;

load('../data/20131114_galdeletions/output/plates_hists_EMD.mat');

data = struct2cell(plates_hists);

d = [0,0;8,12;0,0;0,12;16,12;16,0];

plates = {'Gal1' 'Gal2' 'Gal3' 'Gal4' 'Gal80'};

gal_final = [0 2.^(-8:2)]
glc_final =[0 2.^(-7:-1)];

gal_m = [2.^(-10:3)];glc_m =[2.^(-9:0)];
x_tick = log2(gal_m);
y_tick = log2(glc_m);

%% Plot heatmap MT
for i = [1:5]
    [E_area{i},E_prec{i},E_mean{i}] = Plates2mat({plates{i}},data(i),plates_hists,d,map,th_const,off_peak);
    [D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
    
    figure(i)
    h = pcolor(log2(gal_m),log2(glc_m),M_area{i});hold on;
    set(h,'edgecolor','none');
    
    %     [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,cutoff,gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
    %     plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1]*0,'color',[0 0 0],'markersize',6);hold on
    %     plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
    %
    axis square
    set(gca,'xtick',x_tick([1  3 5:2:end])+0.5,'xticklabel',['No',  mat2cell(-8:2:2)]);
    set(gca,'ytick',y_tick([1  3 5:2:end])+0.5,'yticklabel',['No',mat2cell(-7:2:-1)]);
    colormap(cmap)
    title(['mut',plates{i}]);
    
    Set_fig_YS(figure(i),17,12,12)
    %     print('-dtiff','-r1000',[figure_folder,'mt',plates{i}]);
    
end

%% Plot heatmap WT
for i = [1:5]
    [E_area{i},E_prec{i},E_mean{i}] = Plates2matMch({plates{i}},data(i),plates_hists,d,map,th_const,off_peak);
    [D_area_wt{i},M_area_wt{i}] = ParseHeatmapMat(E_area{i});
    
    figure(i+5)
    h = pcolor(log2(gal_m),log2(glc_m),M_area_wt{i});hold on;
    set(h,'edgecolor','none');
    
    %     [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,cutoff,gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
    %     plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1]*0,'color',[0 0 0],'markersize',6);hold on
    %     plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
    %
    axis square
    set(gca,'xtick',x_tick([1  3 5:2:end])+0.5,'xticklabel',['No',  mat2cell(-8:2:2)]);
    set(gca,'ytick',y_tick([1  3 5:2:end])+0.5,'yticklabel',['No',mat2cell(-7:2:-1)]);
    colormap(cmap)
    title(['wt',plates{i}]);
    
    Set_fig_YS(figure(i+5),17,12,12)
    %     print('-dtiff','-r1000',[figure_folder,'wt',plates{i}]);
    
end

%% plot S-plots
close  all

color_vec = [0 0 0;1 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
gal= gal_final(2:end);
glc = glc_final(2:end);

for i = 1:length(glc)
    for j = 1:length(gal)
        
        rat(i,j) = gal(j)/glc(i);
        
    end
end

figure(6)

j=1;
for i = [1:5]
    figure(i)
    plot(rat,D_area_wt{i},'o','color',color_vec(1,:),'markersize',2,'markerfacecolor',color_vec(1,:));hold on;
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D_area_wt{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(1,:),'linewidth',4);
    
    plot(rat,D_area{i},'o','color',color_vec(2,:),'markersize',2,'markerfacecolor',color_vec(2,:));hold on;
    
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D_area{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(2,:),'linewidth',4);
    
    
    j=j+1;
    box off;xlim([10^-3,10^3.4]);ylim([0 1]);
    Set_fig_YS(figure(i),18,12,12);
end


%% plot decision fronts
close all
fit_cuttoff = [2^-2 2^-8];
cutoff=0.2;
j=1;
for i = [1:5]
    figure(1)
    
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{i},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
    plot(log2(x),log2(y),'o','markerfacecolor',color_vec(j,:),'color',color_vec(j,:),'markersize',2);hold on
    plot(log2(gal),s(log2(gal)),'color',color_vec(j,:),'linewidth',3); hold on;
    
    ylim([-6 0]);xlim([-9 2]);axis square
    j=j+1;
end

Set_fig_YS(figure(1),17,12,12)

