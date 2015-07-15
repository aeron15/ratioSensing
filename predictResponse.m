function predictResponse()

%% Predicted response from multiplication of the single response. Figure S7.

number = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

for i = 1:8
    for j = 1:12
        map{i,j} = [well{i} number{j}];
    end
end

th_const = 3;
off_peak = 2;


cmap=cbrewer('seq', 'YlOrRd', 25);
%figure_folder='/Users/rae10/Dropbox/Galactose_Pathway/gal_paper/Matlab/Matlab_previous_code/MatlabForRenan'
%figure_folder = '../../../figures/Matlab Source';
figure_folder = '.';

%% s288c replicate 2


gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-9:0.5:0]];


gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-10:0.5:0.5];

% load  Replicate 2
plates = {'Plate_A' 'Plate_B' 'Plate_C' 'Plate_D' 'Plate_E' 'Plate_F'};
d = [0,0;0,12;8,0;8,12;16,0;16,12];
%load('C:\Users\RE151\Dropbox\gal_paper\Data\s288c\20140217_4283\output\plates_hists_EMD')
load('../data/s288c/20140217_4283/output/plates_hists_EMD');
data = struct2cell(plates_hists);


[E_area{1},E_prec{1},E_mean{1}] = Plates2mat(plates,data,plates_hists,d,map,th_const,off_peak);

% remove low and high gloucse rows:

E_area{1} = E_area{1}(1:end-4,:); %remove NaN from high glc
E_prec{1} = E_prec{1}(1:end-4,:); %remove NaN from high glc
E_mean{1} = E_mean{1}(1:end-4,:); %remove NaN from high glc

i = 1;
[D_area{i},M_area{i}] = ParseHeatmapMat(E_area{i});
[D_prec{i},M_prec{i}] = ParseHeatmapMat(E_prec{i});
[D_mean{i},M_mean{i}] = ParseHeatmapMat(E_mean{i});
M_area{i}(7,3) = mean([M_area{i}(6,3),M_area{i}(8,3)]);


%% single response

gal_r_area = E_area{1}(1,:);
glc_r_area= E_area{1}(:,end-2); % 2% GAL

figure(4)
subplot(2,1,1)
plot(log2(gal_final),gal_r_area,'o-k');hold on;
xlim([-9.5 2.5]);set(gca,'xtick',[-9:3:2,2]);
ylim([0 1])
subplot(2,1,2)
plot(log2(glc_final),glc_r_area,'o-k');hold on;
xlim([-9.5 0.5]);set(gca,'xtick',[-9:3:0]);
ylim([0 1])

Set_fig_RE(figure(4),17,12,12)
%%

for i = 1:length(gal_r_area)
    for j = 1:length(glc_r_area)
        R(j,i) = gal_r_area(i)*glc_r_area(j);
    end
end
[D_R,M_R] = ParseHeatmapMat(R);

figure(1)
h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-9:1:0)])

colormap(cmap);

figure(2)
h = pcolor(M_area{1});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-9:1:0)])

colormap(cmap)

Set_fig_RE(figure(1),17,12,12)
Set_fig_RE(figure(2),17,12,12)

%% >>>>>>>>>>>>> A10 NEED TO FIND THIS DATA. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
clear R
gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-10:1:1]];

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-12:1:2];

th_const = 2.7;
i=1;
glc = [-Inf   -10    -9    -8    -7    -6    -5    -4    -3    -2    -1     0     1]


load('../data/A10H3/output/plates_hists_EMD_stats');

plates = {'A10a__H3a' 'A10b__H3b' 'A10c__H3c' 'A10d__H3d'};

data = struct2cell(plates_hists);


d = [8,0;8,12;0,0;0,12;16,12;16,0];

data = struct2cell(plates_hists);
i=2;
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

%% single response

fig3=figure(3)
scrsz = get(0,'ScreenSize');
set(fig3,'Position',[1 scrsz(4)*0.45 scrsz(3)*0.45 scrsz(4)])
gal_r_area = E_area{2}(1,:);
glc_r_area= E_area{2}(:,end-2);

figure(3)
subplot(2,1,1)
plot(log2(gal_final),gal_r_area,'o-k','linewidth',2);hold on;
xlim([-9.5 2.5]);set(gca,'xtick',[-9:3:2,2]);
ylim([0 1])


subplot(2,1,2)
plot(log2(glc_final),glc_r_area,'o-k','linewidth',2);hold on;
xlim([-10 1.5]);set(gca,'xtick',[-10:3:1,1]);
ylim([0 1])

Set_fig_RE(figure(3),18,18,18)

filename=['Single_response_A10'];
export_fig(filename, '-pdf','-transparent','-nocrop');
%%
for i = 1:length(gal_r_area)
    for j = 1:length(glc_r_area)
        R(j,i) = gal_r_area(i)*glc_r_area(j);
    end
end
[D_R,M_R] = ParseHeatmapMat(R);

%%
fig1=figure(1)
scrsz = get(0,'ScreenSize');
set(fig1,'Position',[1 scrsz(4)*0.45 scrsz(3)*0.45 scrsz(4)])

h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])

colormap(cmap);

Set_fig_RE(figure(1),18,18,18)


filename=['Predicted_response_A10'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%%


fig2=figure(2)
scrsz = get(0,'ScreenSize');
set(fig2,'Position',[1 scrsz(4)*0.45 scrsz(3)*0.45 scrsz(4)])

h = pcolor(M_area{2});

set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])

colormap(cmap)

Set_fig_RE(figure(2),18,18,18)

filename=['Measured_response_A10'];
export_fig(filename, '-pdf','-transparent','-nocrop');


%% Deplition estimation

close all

clear R t g gg r_rect r_ratio
cutoff_gal = -7;
cutoff_glc = -7;
od0 = 0.003;
d = 0.02;

factor = [1 5 10];


gal_final = [0 2.^[-9:0.5:2]];
glc_final = [0 2.^[-10:1:1]];

gal_m = 2.^[-10:0.5:2.5];
glc_m = 2.^[-12:1:2];

g0 = log(2)/(90/60);

% single
n = 3;
gal_r = 1./(1 + (2^cutoff_gal./gal_final).^n); 
glc_r = 1./(1 + (glc_final/2^cutoff_glc).^n);

gal_r = zeros(size(gal_final))
glc_r = zeros(size(glc_final))

gal_r(gal_final>=2^cutoff_gal) = 1;

glc_r(glc_final<=2^cutoff_glc) = 1;


plot(log2(gal_final),gal_r);
figure
plot(log2(glc_final),glc_r);
%multiplicative
for i = 1:length(gal_final)
    for j = 1:length(glc_final)
        r_rect(j,i) = gal_r(i)*glc_r(j);% [g_glc(j)];
    end
end
[D_R,M_R] = ParseHeatmapMat(r_rect);
figure(1)
h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No',' ', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No',' ', mat2cell(-10:2:1)]);
colormap(cmap);
Set_fig_RE(figure(1),17,12,12)

% ratio 

for i = 1:length(gal_final)
    for j = 1:length(glc_final)
        
        ratio_sep = 2^cutoff_glc*(-1 + gal_final(i)/ 2^cutoff_gal);
        if glc_final(j)>ratio_sep
            r_ratio(j,i) = 0;
        else
            r_ratio(j,i) = 1;
        end
    end
end
[D_R,M_R] = ParseHeatmapMat(r_ratio);
figure(2)
h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1 2 3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])
colormap(cmap);
Set_fig_RE(figure(2),17,12,12)

% growth rates
% g_glc = [0 mat2gray(log2(glc_final(2:end)))]*g0;
% g_gal = [0 mat2gray(log2(gal_final(2:end)))]*g0;

g_glc = g0 * glc_final.^2./(glc_final.^2 + (2^-6)^2);
g_gal = g0 * gal_final.^2./(gal_final.^2 + (2^-6)^2);
% 
% g_glc = [0 mat2gray((glc_final(2:end)))]*g0;
% g_gal = [0 mat2gray((gal_final(2:end)))]*g0;

for i = 1:length(gal_final)
    for j = 1:length(glc_final)
%         g(j,i) = max([g_glc(j),g_gal(i)]);
%         g(j,i) = mean([g_glc(j),g_gal(i)]);
        g(j,i) = g0% [g_glc(j)];
    end
end
% g = g0 * glc_final.^2./(glc_final.^2 + (2^-6)^2);
%% rect deplition
od = od0*[0 1 5 50];

close all
for k =1:length(od)
    clear t R
    for i = 1:length(gal_final)
        for j = 1:length(glc_final)
            
            
            c0 = glc_final(j);
            ct = 2^cutoff_glc;
            if c0<ct
                t(j,i) = 0;
            else
                t(j,i) = 1/g(j,i) * log(  (c0 - ct +d*od(k)) /(d*od(k)) );
            end
            
            if t(j,i)<8 && gal_final(i)>=2^cutoff_gal
                R(j,i) = min((8-t(j,i))/6,1)*r_rect(1,i);
            else
                R(j,i) = 0;
            end
        end
    end
    [D_R,M_R] = ParseHeatmapMat(R);
    figure(k)
    h = pcolor(M_R);
    set(h,'edgecolor','none');
    axis square
    set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])
    colormap(cmap);
    Set_fig_RE(figure(k),17,12,12)
    
    
    [D_R,M_R] = ParseHeatmapMat(t);
    figure(10*k)
    h = pcolor(M_R);
    set(h,'edgecolor','none');
    axis square
    set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
    set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])
    colormap(cmap);
    Set_fig_RE(figure(k),17,12,12)
end
%%
% ratio deplition
close all
for k =1:length(od)

for i = 1:length(gal_final)
    for j = 1:length(glc_final)
        
        
        c0 = glc_final(j);
        ratio_sep = 2^cutoff_glc*(-1 + gal_final(i)/ 2^cutoff_gal);
        ct = ratio_sep;
        if c0<ct
            t(j,i) = 0;
        else
        t(j,i) = 1/g(j,i) * log(  (c0 - ct + d*od(k)) /(d*od(k)) );
        end
        
        if t(j,i)<8 && gal_final(i)>=2^cutoff_gal
            R(j,i) = min((8-t(j,i))/6,1)*r_ratio(1,i);
        else
            R(j,i) = 0;
        end
    end
end
[D_R,M_R] = ParseHeatmapMat(R);
figure(100*k)
h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])
colormap(cmap);
Set_fig_RE(figure(100*k),17,12,12)


[D_R,M_R] = ParseHeatmapMat(t);
figure(1000*k)
h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-10:2:1)])
colormap(cmap);
Set_fig_RE(figure(1000*k),17,12,12)
end