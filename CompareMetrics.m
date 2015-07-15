function CompareMetrics()
%% compare metrics

number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end

th_const = 3;
off_peak = 2;


gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-9:0.5:0]];

cmap=cbrewer('seq', 'YlOrRd', 25);
%figure_folder = 'C:\Users\ys151\Dropbox\gal_paper\figures\Matlab Source\';
figure_folder = '	';

%% load  Replicate 2
plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
load('../data/s288c/20140217_4283/output/plates_hists_EMD')
data = struct2cell(plates_hists);

[E_area{2},E_prec{2},E_mean{2}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

E_area{2} = E_area{2}(1:end-4,:); %remove NaN from high glc
E_prec{2} = E_prec{2}(1:end-4,:); %remove NaN from high glc
E_mean{2} = E_mean{2}(1:end-4,:); %remove NaN from high glc

%
i = 2;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_prec{i},M_prec{i}] = ParseHeatmapMat(E_prec{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});

D{1} = D_area{2};
D{2} = D_prec{2};

D{3} = (D_mean{2});
D{3}(find(isinf(D{3})))=nan;
D{3}=D{3}-min(min(D{3}));
D{3} = D{3}/(max(max(D{3})));

M{1} = M_area{2};
M{2} = M_prec{2};
M{3} = (M_mean{2});

th_const_vec = [2.5 3 3.5];

for k = 1:length(th_const_vec)
    
    [E_area{2},E_prec{2},E_mean{2}] = Plates2mat(plates,data,plates_hists,d,map,th_const_vec(k),off_peak);
    
    % remove low and high gloucse rows:
    
    E_area{2} = E_area{2}(1:end-4,:); %remove NaN from high glc
    E_prec{2} = E_prec{2}(1:end-4,:); %remove NaN from high glc
    E_mean{2} = E_mean{2}(1:end-4,:); %remove NaN from high glc
    
    %
    i = 2;
    [D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
    [D_prec{i},M_prec{i}] = ParseHeatmapMat(E_prec{i});
    [D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});
    
    %Store the next 3 different thresholds (2.5, 3 and 3.5 for the
    %percentage metric in the positions 4,5 and 6 but this could be
    %different
    
    D{k+3} = D_prec{2};
    M{k+3} = M_prec{2};
    
    
end

%%
i
figure(1)
h = pcolor(M_area{i});
title('Area')
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

figure(2)
h = pcolor(M_prec{i});
title('Threshold')
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

figure(3)
h = pcolor(M_mean{i});
title('Mean')
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

%% Pull out specific histogram of high and low induction

%Facs

% from replicate 2 (-6.5,-8) (-6.5,-2)
p = 1; i=1;j = 1;


plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
%load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\20140217_4283\output\plates_hists_EMD')
load('../data/s288c/20140217_4283/output/plates_hists_EMD')

data = struct2cell(plates_hists);

plate_name = 'Plate_A';well_name = 'A01';
[x,y] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_B';well_name = 'A12';
[x2,y2] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

figure(1)
plot(x,y,'r');hold on;plot(x2,y2,'k')

%% Make the S-plots for different mesures
close all
gal= gal_final(2:end);
glc = glc_final(2:end);

clear rat
for i = 1:length(glc)
    for j = 1:length(gal)
        
        rat(i,j) = gal(j)/glc(i);
        
    end
end

color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];

figure(6)
for i = [1 2 3]
    %     figure(i)
    plot(rat,D{i},'o','color',color_vec(i,:),'markersize',2,'markerfacecolor',color_vec(i,:));hold on;
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    %     plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(i,:),'linewidth',4)
end
box off;xlim([10^-3,10^3.4]);ylim([0 1]);
Set_fig_RE(figure(6),18,12,12);
%% Plot Contour for different mesures
close all

fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-10:0.5:0.5];

x_tick = log2(gal_m);
y_tick = log2(glc_m);


th_vec = [2.5 3 3.5];
prec_vec = [0.2 0.35 0.5 0.65 0.8 0.9]
for i = [1 2 3 4 5 6]
    
    for j = 1:length(prec_vec)
        M{i}(7,3) = mean([M{i}(6,3),M{i}(8,3)]);
        
            M_area{i}(1,1) = 1;
        figure(i)
        h = pcolor(log2(gal_m),log2(glc_m),M{i});hold on;
        set(h,'edgecolor','none');
        
        
        axis square
        set(gca,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
        set(gca,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
            h = colorbar('location','northoutside');
        colormap(cmap)
        
        figure(10*i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,prec_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        plot(log2(x),log2(y),'o','markerfacecolor',color_vec(j,:),'color',color_vec(j,:),'markersize',4);hold on
            plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
        axis square
        xlim([-10 2.5]);ylim([-10 0.5]);
        
            set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
        Set_fig_RE(figure(i),17,12,12)
        Set_fig_RE(figure(10*i),17,12,12)
        
    end
        print('-dtiff','-r1000',[figure_folder,'measures',num2str(i)]);
    
end
%% load A10
gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-10:1:1]];
th_const = 2.7;
i=1;
glc = [-Inf   -10    -9    -8    -7    -6    -5    -4    -3    -2    -1     0     1]

load('../data/A10H3/output/plates_hists_EMD_stats')

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
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-10:2:1)])
axis square
colorbar
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

M{1} = M_area{1};
M{2} = M_prec{1};
M{3} = M_mean{1};

D{1} = D_area{1};
D{2} = D_prec{1};
D{3} = (D_mean{1});D{3}(find(isinf(D{3})))=nan;D{3}=D{3}-min(min(D{3}));D{3} = D{3}/(max(max(D{3})));



%% different threshold


th_const_vec = [2.5 2.75 3 3.5];

i=4
for k = 1:length(th_const_vec)
    
    [E_area{i},E_prec{i},E_mean{i}] = Plates2mat(plates,data,plates_hists,d,map,th_const_vec(k),off_peak);
    
    
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
    
    D{k+3} = D_prec{i};
    M{k+3} = M_prec{i};
    
end


%% Plot Contour for different mesures
fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-12:1:2];

x_tick = log2(gal_m);
y_tick = log2(glc_m);


prec_vec = [0.2 0.35 0.5 0.65 0.8 0.95]
for i = [1 2 3]
    for j = 1:length(prec_vec)
        
        figure(i)
        h = pcolor(log2(gal_m),log2(glc_m),M{i});hold on;
        set(h,'edgecolor','none');
        colormap(cmap)
        axis square
        set(gca,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
        set(gca,'ytick',y_tick([1  3 5:2:13])+0.5,'yticklabel',['No', mat2cell(-10:2:1)]);
        
        figure(10*i)
        [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,prec_vec(j),gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
        plot(log2(x),log2(y),'o','markerfacecolor',color_vec(j,:),'color',color_vec(j,:),'markersize',4);hold on
        %     plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
        axis square
        xlim([-10 2.5]);ylim([-10 0.5]);
        
        
        Set_fig_RE(figure(i),17,12,12)
        Set_fig_RE(figure(10*i),17,12,12)
        
    end
    %     print('-dtiff','-r1000',[figure_folder,'measures',num2str(i)]);
    
end

%%

%% Pull out specific histogram of high and low induction for A10

%Facs

% from replicate 2 (-6.5,-8) (-6.5,-2)
p = 1; i=1;j = 1;


plates = {'A10a__H3a' 'A10b__H3b' 'A10c__H3c' 'A10d__H3d'};

d = [0,0;0,12;8,0;8,12;16,0;16,12];
load('../data/A10H3/output/plates_hists_EMD_stats')

data = struct2cell(plates_hists);

plate_name = 'A10a__H3a';well_name = 'A01';
[x,y] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'A10d__H3d';well_name = 'A12';
[x2,y2] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

figure(1)
plot(x,y,'r');hold on;plot(x2,y2,'k')
%%
%% Make the S-plots for different mesures
close all
gal= gal_final(2:end);
glc = glc_final(2:end);

clear rat
for i = 1:length(glc)
    for j = 1:length(gal)
        
        rat(i,j) = gal(j)/glc(i);
        
    end
end

color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];

figure(6)
for i = [1 2 3]
    %     figure(i)
    plot(rat,D{i},'o','color',color_vec(i,:),'markersize',4,'markerfacecolor',color_vec(i,:));hold on;
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    %     plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(i,:),'linewidth',4)
end
box off;xlim([10^-3,10^3.4]);ylim([0 1]);
Set_fig_RE(figure(6),18,12,12);

end