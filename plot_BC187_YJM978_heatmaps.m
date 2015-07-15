%% load A10
% created by RE 20141006


number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end

off_peak = 2;


cmap=cbrewer('seq', 'YlOrRd', 25);
figure_folder = '	';


gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-10:1:1]];
th_const = 2.7;
i=1;
glc = [-Inf   -10    -9    -8    -7    -6    -5    -4    -3    -2    -1     0     1]

load('../../../Data/A10H3/plates_hists_EMD_stats')

plates = {'A10a__H3a' 'A10b__H3b' 'A10c__H3c' 'A10d__H3d'};

data = struct2cell(plates_hists);


d = [8,0;8,12;0,0;0,12;16,12;16,0];

data = struct2cell(plates_hists);

[E_area{i},E_prec{i},E_mean{i}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% add nan for the nonexisiting glc low levels

temp = nan*ones(13,24);
temp(1,:) = E_area{i}(1,:);
temp2(1,:) = E_prec{i}(1,:);
temp3(1,:) = E_mean{i}(1,:);

for j = 1:length(glc_final)
    ind = find(log2(glc_final(j))==glc);
    if ~isempty(ind)
        temp(j,:) = E_area{i}(ind,:);
        temp2(j,:) = E_prec{i}(ind,:);
        temp3(j,:) = E_mean{i}(ind,:);
    end
end

E_area{i}= temp;
E_prec{i}= temp2;
E_mean{i}= temp3;

[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_prec{i},M_prec{i}] = ParseHeatmapMat(E_prec{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(i)
%h = pcolor(M_area{i});
%h = pcolor(M_area{i}(:,([1:23 26])));
h = pcolor(M_mean{i}(:,([1:23 26])));

set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-10:2:1)])
axis square
colorbar('northoutside')
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

M{1} = M_area{1};
M{2} = M_prec{1};
M{3} = M_mean{1};

D{1} = D_area{1};
D{2} = D_prec{1};
D{3} = (D_mean{1});D{3}(find(isinf(D{3})))=nan;D{3}=D{3}-min(min(D{3}));D{3} = D{3}/(max(max(D{3})));

filename=['BC187_heatmap_mean'];
export_fig(filename, '-pdf','-transparent','-nocrop');
close all;

%% load H3


data = struct2cell(plates_hists);


d = [8,0;8,12;0,0;0,12;16,12;16,0];

data = struct2cell(plates_hists);

[E_area{i},E_prec{i},E_mean{i}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);

% add nan for the nonexisiting glc low levels

temp = nan*ones(13,24);
temp(1,:) = E_area{i}(1,:);
temp2(1,:) = E_prec{i}(1,:);
temp3(1,:) = E_mean{i}(1,:);

for j = 1:length(glc_final)
    ind = find(log2(glc_final(j))==glc);
    if ~isempty(ind)
        temp(j,:) = E_area{i}(ind,:);
        temp2(j,:) = E_prec{i}(ind,:);
        temp3(j,:) = E_mean{i}(ind,:);
    end
end

E_area{i}= temp;
E_prec{i}= temp2;
E_mean{i}= temp3;

[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_prec{i},M_prec{i}] = ParseHeatmapMat(E_prec{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

figure(i)
%h = pcolor(M_area{i});
% h = pcolor(M_area{i}(:,([1:23 26])));
h = pcolor(M_mean{i}(:,([1:23 26])));

set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-10:2:1)])
axis square
colorbar('northoutside')
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

M{1} = M_area{1};
M{2} = M_prec{1};
M{3} = M_mean{1};

D{1} = D_area{1};
D{2} = D_prec{1};
D{3} = (D_mean{1});D{3}(find(isinf(D{3})))=nan;D{3}=D{3}-min(min(D{3}));D{3} = D{3}/(max(max(D{3})));


filename=['YJM978_heatmap_mean'];
export_fig(filename, '-pdf','-transparent','-nocrop');
close all;

