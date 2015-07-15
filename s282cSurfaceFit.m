function s282cSurfaceFit

%% Collect and combine all the s282c replictaes into one format

% M matrix for ploting including zeros in a sperate row
% D matrix for fitting purposes

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

cmap='jet';%cbrewer('seq', 'YlOrRd', 25);
%figure_folder = 'C:\Users\ys151\Dropbox\gal_paper\figures\Matlab Source\';
figure_folder = '.';
% cmap = buildcmap_YS([254,237,222;253,190,133;253,141,60;217,71,1]/255);

%% Replicate 1
%load('C:\Users\yoni\Dropbox\gal_paper\Data\Gal2FullD\20140331_large_DG_gal2del\output\plates_hists')
load('../data/20140331_large_DG_gal2del/output/plates_hists')

gal = [0 2.^[-9:0.5:2]];
glc = [0 2.^[-10.5:0.5:0]];

d = [0,0;0,12;8,0;8,12;16,0;16,12];

data = struct2cell(plates_hists);
plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};


[E_area{1},E_prec{1},E_mean{1}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:
E_area{1} = E_area{1}(1:end-1,:); %remove NaN from high glc
E_area{1} = E_area{1}([1 5:end],:);% remove low glc

E_prec{1} = E_prec{1}(1:end-1,:); %remove NaN from high glc
E_prec{1} = E_prec{1}([1 5:end],:);% remove low glc

E_mean{1} = E_mean{1}(1:end-1,:); %remove NaN from high glc
E_mean{1} = E_mean{1}([1 5:end],:);% remove low glc

%
i = 1;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar 
colormap(cmap)

Set_fig_RE(figure(1),12,12,12)
% figure(6)
% plot(x_off_bfp,y_off_bfp);hold on
%% Replicate 2
plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
%load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\20140217_4283\output\plates_hists_EMD')
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

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)
% figure(6)
% plot(x_off_bfp,y_off_bfp);hold on
%% Replicate 3

plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};

d = [0,0;0,12;8,0;8,12;16,12;16,0];

load('../data/s288c/20140217_YM4277_mig1_del_Gal4/output/plates_hists_EMD')

data = struct2cell(plates_hists);

[E_area{3},E_prec{3},E_mean{3}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

E_area{3} = E_area{3}(1:end-4,:); %remove NaN from high glc


i = 3;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});


figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),17,12,12)

%% Replicate 4

load('../data/s288c/20131019_2D_LG/output/plates_hists_EMD')

plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D'};

data = struct2cell(plates_hists);

gal = [0 2.^[-9:0.5:2]];                                    0 
glu = [0 2.^[-7:0.5:0]];

d = [0,0;0,12;8,0;8,12;16,12;16,0];

data = struct2cell(plates_hists);

[E_area{4},E_prec{4},E_mean{4}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% add nan for the nonexisiting glc low levels

temp = nan*ones(size(E_area{3}));
temp(1,:) = E_area{4}(1,:);
temp(6:end,:) = E_area{4}(end-14:end,:);
E_area{4} = temp;

i = 4;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),17,12,12)

%% A10 and H3

%% A10

glc = [-Inf   -10    -9    -8    -7    -6    -5    -4    -3    -2    -1     0     1]

load('../data/A10H3/output/plates_hists_EMD_stats')

plates = {'A10a__H3a' 'A10b__H3b' 'A10c__H3c' 'A10d__H3d'};

data = struct2cell(plates_hists);


d = [8,0;8,12;0,0;0,12;16,12;16,0];

data = struct2cell(plates_hists);

[E_area{5},E_prec{5},E_mean{5}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% add nan for the nonexisiting glc low levels

temp = nan*ones(20,24);
temp(1,:) = E_area{5}(1,:);
for i = 1:length(glc_final)
    ind = find(log2(glc_final(i))==glc);
    if ~isempty(ind)
    temp(i,:) = E_area{5}(ind,:);
    end
end

E_area{5} = temp;

i = 5;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

%% H3

plates = {'A10a__H3a' 'A10b__H3b' 'A10c__H3c' 'A10d__H3d'};

data = struct2cell(plates_hists);

d = [8,0;8,12;0,0;0,12;16,12;16,0];

data = struct2cell(plates_hists);

[E_area{6},E_prec{6},E_mean{6}] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak);

% add nan for the nonexisiting glc low levels

temp = nan*ones(20,24);
temp(1,:) = E_area{6}(1,:);
for i = 1:length(glc_final)
    ind = find(log2(glc_final(i))==glc);
    if ~isempty(ind)
    temp(i,:) = E_area{6}(ind,:);
    end
end

E_area{6} = temp;

i = 6;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});

figure(i)
h = pcolor(M_area{i});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar
colormap(cmap)

Set_fig_RE(figure(i),12,12,12)

%% Median of s282c

for i = 1:4
    MM(:,:,i) = M_area{i};
end

med = nanmedian(MM,3);
figure(5)
h = pcolor(med);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-9:1:0)])

colorbar('location','northoutside')
colormap(cmap)

Set_fig_RE(figure(5),17,12,12)

%% Make the S-plots. For replcate 2 Data/s288c/20140217_4283/output/plates_hists_EMD

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
for i = [2]
%     figure(i)
    plot(rat,D_area{i},'o','color',color_vec(i,:),'markersize',4,'markerfacecolor',color_vec(i,:));hold on;
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D_area{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    %n(i)=s.n;
    %c(i)=s.c;
    plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(i,:),'linewidth',4)
end
box off;xlim([10^-3,10^3.4]);ylim([0 1]);
% N = mean(n(2:end));dN=std(n(2:end));
% C = mean( c(2:end).^(1./n(2:end)));dC = std( c(2:end).^(1./n(2:end)))
ylabel('Fraction of cells inducing');
xlabel('log_{10} (Galactose/Glucose)')
Set_fig_RE(figure(6),18,20,20);


filename=['S_plot_S288C'];
% export_fig(filename, '-pdf','-nocrop');

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
D{1} = D_area{2};
D{2} = D_prec{2};
D{3} = (D_mean{2});D{3}(find(isinf(D{3})))=nan;D{3}=D{3}-min(min(D{3}));D{3} = D{3}/(max(max(D{3})));
figure(6)
for i = [1 2 3]
%     figure(i)
    plot(rat,D{i},'o','color',color_vec(i,:),'markersize',2,'markerfacecolor',color_vec(i,:));hold on;
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(i,:),'linewidth',4)
end
box off;xlim([10^-3,10^3.4]);ylim([0 1]);
Set_fig_RE(figure(6),18,12,12);
%% Plot Contour for different mesures
fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-10:0.5:0.5];

x_tick = log2(gal_m);
y_tick = log2(glc_m);

M{1} = M_area{2};
M{2} = M_prec{2};
M{3} = (M_mean{2});

for i = [1 2 3]
 figure(i)

     M{i}(7,3) = mean([M{i}(6,3),M{i}(8,3)]);

%     M_area{i}(1,1) = 1;
    h = pcolor(log2(gal_m),log2(glc_m),M{i});hold on;
    set(h,'edgecolor','none');
    
%     [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D{i},1,0,cutoff,gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
%     plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1]*0,'color',[0 0 0],'markersize',6);hold on
%     plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
%     
    axis square
    set(gca,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(gca,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
%     h = colorbar('location','northoutside');
        colormap(cmap)

%     set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
    Set_fig_RE(figure(i),17,12,12)
 
end

% print('-dtiff','-r1000',[figure_folder,'1B']);
%% Plot Contour 
fit_cuttoff = [2^-2 2^-9];
mid_value = 2^-4;
cutoff=0.2;
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-10:0.5:0.5];

x_tick = log2(gal_m);
y_tick = log2(glc_m)

for i = [2 3 4]
 figure(i)
 if i == 2
     M_area{2}(7,3) = mean([M_area{2}(6,3),M_area{2}(8,3)]);
 end
%     M_area{i}(1,1) = 1;
    h = pcolor(log2(gal_m),log2(glc_m),M_area{i});hold on;
    set(h,'edgecolor','none');
    
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D_area{i},1,0,cutoff,gal_final(2:end),glc_final(2:end),fit_cuttoff,mid_value);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1]*0,'color',[0 0 0],'markersize',6);hold on
    plot(log2(x),s(log2(x)),'color',[0 0 0],'linewidth',2); hold on;
    
    axis square
    set(gca,'xtick',x_tick([1  3 5:2:25])+0.25,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(gca,'ytick',y_tick([1  3 5:2:22])+0.25,'yticklabel',['No', mat2cell(-9:1:0)]);
%     h = colorbar('location','northoutside');
        colormap(cmap)

%     set(h, 'XTick',[0:0.2:1],'XTicklabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'});
    Set_fig_RE(figure(i),17,12,12)
 
end

% print('-dtiff','-r1000',[figure_folder,'1B']);
%% plot decision fronts
close all
fit_cuttoff = [2^-2 2^-8];
cutoff=0.2;
mid_value = 2^-4;

for i = [2]
 figure(1)
 if i>=5
     D = D_area{i}(1:2:19,:);
     gal = gal_final(2:end);
     glc = glc_final(2:2:end);
 else
     D = D_area{i};
     gal = gal_final(2:end);
     glc = glc_final(2:end);
 end
%     M_area{i}(1,1) = 1;
%     h = pcolor(log2(gal_m),log2(glc_m),M_area{i});hold on;
%     set(h,'edgecolor','none');
    
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(D,1,0,cutoff,gal,glc,fit_cuttoff,mid_value);
    plot(log2(x),log2(y),'o','markerfacecolor',color_vec(i,:),'color',color_vec(i,:),'markersize',2);hold on
    plot(log2(gal),s(log2(gal)),'color',color_vec(i,:),'linewidth',6); hold on;
    
    ylim([-9 0]);xlim([-9 2]);axis square
    Set_fig_RE(figure(1),17,12,12)

end

%% Pure ratio
cmap=cbrewer('seq', 'YlOrRd', 25);


close all
figure(10)
subplot(1,2,1)
a = logspace(-2,3,100);
b = logspace(-2,3,100);
[A,B] = meshgrid(a,b);

    ff = (1./(1 + (B./A).^2));
    h= contourf(a,b,ff,5,'edgecolor','none');hold on;
       c = [0.01 1 10 1000]
    for i = 1:length(c)
    plot(a,c(i)*1./a,'k'); 
    end
    set(gca,'xscale','log','yscale','log');

    axis('square');view([0 90]);
%     xlim([10^2 10^7]);ylim([10^2 10^7])
%     xlabel('A');ylabel('R');
    colormap(cmap);caxis = [0 1];
    box off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    
    Set_fig_RE(figure(10),12,18,18)
print('-dtiff','-r1000',[figure_folder,'1D1']);


    for i = 1:length(b)
        for j = 1:length(a)
            
            rat(i,j) = a(j)/b(i);
            
        end
    end
    subplot(1,2,2)
%
        plot(rat,ff,'color','k','linewidth',2);hold on;
    set(gca,'xscale','log');
%     [XOut, YOut, ZOut] = prepareSurfaceData(a, b, ff);
%     s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
%     plot(rat,s.c./(s.c + (rat.^-s.n)),'color','k','linewidth',2)
         box off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
        xlim([10^-4,10^4]);ylim([-0.01 1.01])
axis('square');

%% threshold 
cmap=cbrewer('seq', 'YlOrRd', 25);


close all
figure(10)
subplot(1,2,1)
m_min = -2;
m_max = 3
a = logspace(m_min,m_max,100);
b = logspace(m_min,m_max,100);
[A,B] = meshgrid(a,b);

    ff = (1./(1 + (1./A).^2)).*(1./(1 + (B./1).^2));
    
    h= contourf(a,b,ff,5,'edgecolor','none');hold on;
    c = [0.01 1 10 1000]
    for i = 1:length(c)
    plot(a,c(i)*1./a,'k'); 
    end
    set(gca,'xscale','log','yscale','log');
    
    axis('square');view([0 90]);
%     xlim([10^2 10^7]);ylim([10^2 10^7])
%     xlabel('A');ylabel('R');
    colormap(cmap);caxis = [0 1];
    box off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    Set_fig_RE(figure(10),12,18,18)
print('-dtiff','-r1000',[figure_folder,'1D2']);


    for i = 1:length(b)
        for j = 1:length(a)
            
            rat(i,j) = a(j)/b(i);
            
        end
    end
    subplot(1,2,2)
%
    
    for i = 1:length(c)
        
            x =  logspace(-4,5,100);
%             x =  logspace(log10(c(i)/(max(b)^2)),log10(c(i)/(min(b)^2)),100);

    ff = (1./(1 + (1./(sqrt(c(i))*sqrt(x))).^2)).*(1./(1 + ((sqrt(c(i))./sqrt(x))./1).^2));

        plot(x,ff,'color','k','linewidth',2);hold on;
    end

        set(gca,'xscale','log');
%     [XOut, YOut, ZOut] = prepareSurfaceData(a, b, ff);
%     s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
%     plot(rat,s.c./(s.c + (rat.^-s.n)),'color','k','linewidth',2)
         box off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
        xlim([10^-4,10^5]);ylim([-0.01 1.01])
axis('square');




%% Pull out specific histogram of a well

%Facs

% from replicate 2 (-6.5,-8) (-6.5,-2)
p = 1; i=7;j = 4;


plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
load('../data/s288c/20140217_4283/output/plates_hists_EMD')

data = struct2cell(plates_hists);

plate_name = 'Plate_A';well_name = 'G04';
[x,y] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_A';well_name = 'G10';
[x2,y2] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_C';well_name = 'D10';
[x3,y3] = GetHistogram(plates,data,plates_hists,plate_name,well_name);

plate_name = 'Plate_D';well_name = 'D02';
[x4,y4] = GetHistogram(plates,data,plates_hists,plate_name,well_name);



max_x = max([max(x),max(x2),max(x3),max(x4)]);
min_x = min([min(x),min(x2),min(x3),min(x4)]);
xx = linspace(min_x,max_x,150);

yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
yy(isnan(yy))=0;

yy2 = interp1(x2,y2,xx,'linear');yy2(yy2<0)=0;
yy2(isnan(yy2))=0;

yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
yy(isnan(yy))=0;

yy2 = interp1(x2,y2,xx,'linear');yy2(yy2<0)=0;
yy2(isnan(yy2))=0;

yy3 = interp1(x3,y3,xx,'linear');yy3(yy3<0)=0;
yy3(isnan(yy3))=0;

yy4 = interp1(x4,y4,xx,'linear');yy4(yy4<0)=0;
yy4(isnan(yy4))=0;


figure(10)
plot(xx,yy/sum(yy),'k','linewidth',6);box off;
Set_fig_RE(figure(10),25,18,18)
print('-dpdf','-r1000',[figure_folder,'1B1']);

figure(11)
plot(xx,yy2/sum(yy2),'k','linewidth',6);box off;
Set_fig_RE(figure(11),25,18,18)
print('-dtiff','-r1000',[figure_folder,'1B2']);

figure(12)
plot(xx,yy3/sum(yy3),'k','linewidth',6);box off;
Set_fig_RE(figure(12),25,18,18)
print('-dtiff','-r1000',[figure_folder,'1B3']);

figure(13)
plot(xx,yy4,'k','linewidth',6);box off;
Set_fig_RE(figure(13),25,18,18)

