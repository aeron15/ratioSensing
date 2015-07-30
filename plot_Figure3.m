function plot_Figure3()
%PLOT_FIGURE3
%Deletion of GAL2 does not eliminate ratio sensing. Black and red lines are the decision front
%inspired by figure 4

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

load('../data/20140517_gal2del_large2DG/output/plates_hists')

glc = 2.^([-Inf   -7.5:0.5:0]);

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};

%% GAL2 delete
i=1;

[E_area{i},E_prec{i},E_mean{i}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

temp = nan*ones(20,24);
for j = 1:length(glc_final)
    ind = find(log2(glc_final(j))==log2(glc));
    if ~isempty(ind)
        temp(j,:) = E_area{i}(ind,:);
    end
end
E_area{i} = temp;

[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])
ylim([6 20]);
colorbar
colormap(cmap)

Set_fig_RE(figure(i),17,12,12)

filename='Figure3A_gal2delete';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

%% WT
i=2;

[E_area{i},E_prec{i},E_mean{i}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);

temp = nan*ones(20,24);
for j = 1:length(glc_final)
    ind = find(log2(glc_final(j))==log2(glc));
    if ~isempty(ind)
        temp(j,:) = E_area{i}(ind,:);
    end
end
E_area{i} = temp;


[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])
ylim([6 20]);

colorbar
colormap(cmap)

Set_fig_RE(figure(i),17,12,12)
filename='Figure3A_WT';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

%% Decision Fronts of gal2 delete and WT
gal= gal_final(2:end);
glc = glc_final(2:end);
color_vec = [0 0 0;1 0 0;0 0 1;0 1/2 0];


fit_cuttoff = [2^-1 2^-9];
mid_value = 2^-4;
cutoff=0.2;

figure(4)
i=1;
for j = [1 2]
    [x,y,s,~,~,~,~,~,~] = SmoothHeatMap(D_area{j},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
    %[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{j},1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
    
    plot(log2(x),log2(y),'o','markerfacecolor',color_vec(i,:),'color',color_vec(i,:),'markersize',2);hold on
    plot(log2(x),s(log2(x)),'color',color_vec(i,:),'linewidth',4); hold on;
    xlim([-9 0]);ylim([-9 2]);
    i = i+1;
end

axis square

Set_fig_RE(figure(4),17,12,12)
filename='Figure3A_decision_fronts_gal2_and_WT';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

end
