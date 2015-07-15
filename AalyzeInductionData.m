%% Analyze Deletion

clear all
close all
clc



%Yoni
% A10 = load('C:\Users\ys151\Dropbox\gal_paper\data\20130628\MX_mean_A10');
% H3 = load('C:\Users\ys151\Dropbox\gal_paper\data\20130628\MX_mean_H3');


A10 = load('../../../data/20130628/MX_mean_A10');
H3 = load('../../../data/20130628/MX_mean_H3')


wt = load('mean_matrix_bfp_yfp.mat');

A10 = (A10.Mx); A10(find(A10==0))=2;;A10 = mat2gray(A10);

AA10{1} = A10;


H3 = (H3.Mx); H3(find(H3==0))=2;H3 = mat2gray(H3);



wt = flipud(mat2gray(wt.mean_matrix));

gal_s288c = [0 2.^[-9:0.5:2]];
glu_s288c = [0 2.^[-7:0.5:0]];

gal_s288b_4283other = [0 2.^[-9:0.5:2]];
glu_s288b_4283other = [0 2.^[-9:0.5:0]];



gal = [0 2.^(-9:0.5:2)]
glu =[0 2.^(-10:1)];

gal4 = [0 2.^[-9:0.5:2]];
glu4 = [0 2.^[-11:0.5:0]];

%load('C:\Users\ys151\Dropbox\gal_paper\Data\Mutanta\Heatmap_matrices.mat');
load('../../../Data/Mutanta/Heatmap_matrices.mat');

wtt{1} = (heatmaps_wt.Gal1.mean.Matrix2);
wtt{2} = (heatmaps_wt.Gal2.mean.Matrix2);
wtt{3} = (heatmaps_wt.Gal3.mean.Matrix2);
wtt{4} = (heatmaps_wt.Gal4.mean.Matrix2);
wtt{5} = (heatmaps_wt.Gal80.mean.Matrix2);

gal3 = [0 2.^(-8:2)];
glu3 =[0 2.^(-6:0)];



%% s288c
clear wt mt
%x = load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\mean_matrix_bfp_yfp.mat');
x = load('../../../Data/s288c/mean_matrix_bfp_yfp.mat');

WT{1} = flipud(x.mean_matrix);
%load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\20140217_4283\output\plates_hists_EMD')
load('../../../Data/s288c/20140217_4283/output/plates_hists_EMD')

data = struct2cell(plates_hists);

for k = 1:6
    clear wt_temp mt_temp
    
    
    xx = struct2cell(data{k});
    for i = 1:length(xx)
        wt_temp(i) = xx{i}.bfp_yfp.mean;
        mt_temp(i) = xx{i}.other_yfp.mean;
    end
    if k<=4
        wt(:,:,k) =  (reshape(wt_temp,12,8));
        mt(:,:,k) =  (reshape(mt_temp,12,8));
    else
        wt(:,:,k) =  [reshape(wt_temp,12,4),nan*ones(12,4)];
        mt(:,:,k) =  [reshape(mt_temp,12,4),nan*ones(12,4)];
    end
    
    
end
WT{2} = [ wt(:,:,1)',wt(:,:,2)' ; wt(:,:,3)' ,wt(:,:,4)';wt(:,1:4,5)' ,wt(:,1:4,6)'];


% load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\20140217_YM4277_mig1_del_Gal4\output\plates_hists_EMD')
% data = struct2cell(plates_hists);
% for k = 1:6
%
%
%
%         xx = struct2cell(data{k});
%         for i = 1:length(xx)
%                 wt_temp(i) = xx{i}.bfp_yfp.mean;
%                 mt_temp(i) = xx{i}.other_yfp.mean;
%         end
%
%         wt(:,:,k) =  reshape(wt_temp,12,8);
%         mt(:,:,k) =  reshape(mt_temp,12,8);
%
% end
%
% WT{3} = [ wt(:,:,1)' wt(:,:,2)' ; wt(:,:,3)' wt(:,:,4)';wt(:,:,5)' wt(:,:,6)'];


figure(1)
pcolor(log2(gal_s288c),log2(glu_s288c),WT{1})

figure(2)
pcolor(log2(gal_s288b_4283other),log2(glu_s288b_4283other),WT{2})

%% A10 H3 and s288c

close all
cutoff = 0.2;
cmap =  buildcmap_YS([1 1 1;0.5 0 0]);
colormap(cmap)
subplot(1,3,1);hold on
h=pcolor(log2(gal2),log2(glu2),wt);
contour(log2(gal2),log2(glu2),wt,cutoff,'k');
set(h,'edgecolor','none')
axis square;title('s288c');
xlim([-9 1]);ylim([-7 0]);

subplot(1,3,2);hold on
h=pcolor(log2(gal),log2(glu),A10);contour(log2(gal),log2(glu),A10,cutoff,'k');axis square;
set(h,'edgecolor','none')
;title('A10')
xlim([-9 1]);ylim([-7 0]);

subplot(1,3,3);hold on
h = pcolor(log2(gal),log2(glu),H3);contour(log2(gal),log2(glu),H3,cutoff,'k');axis square;
set(h,'edgecolor','none')
;title('H3')
xlim([-9 1]);ylim([-7 0]);


Set_fig_YS(1,18,18,18)
%%






%% Combine s288c H3 and A10
close all
fit_cuttoff = [2^-1 2^-7];
mid_value = 2^-4;
cutoff=0.1;

i=1;
WT{i}(find(WT{i}==-inf))=nan;
WT{i}(find(WT{i}==0))=nan;
[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT{i},max(max(WT{i})),min(min(WT{i})),cutoff,gal_s288c,glu_s288c,fit_cuttoff,mid_value);
figure(2);
plot(log2(x),log2(y),'ok','markersize',1.5);hold on;
plot(log2(x),s(log2(x)),'k');hold on;

i=2
WT{i}(find(WT{i}==-inf))=nan;
WT{i}(find(WT{i}==0))=nan;

[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT{i},max(max(WT{i})),min(min(WT{i})),cutoff,gal_s288b_4283other,glu_s288b_4283other,fit_cuttoff,mid_value);
figure(2);
plot(log2(x),log2(y),'ok','markersize',1.5);
plot(log2(x),s(log2(x)),'k');hold on;


%%

for i = 1
    
    wtt{i}(find(wtt{i}==-inf))=nan;
    min_wt = min(min(wtt{i}));
    max_wt = max(max(wtt{i}));
    [x,y,s,a(i+3),b(i+3),a_d(i+3),a_u(i+3),b_d(i+3),b_u(i+3)] = SmoothHeatMap(wtt{i},max_wt,min_wt,cutoff,gal3,glu3,fit_cuttoff,mid_value);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 0],'markersize',1.5);
    figure(2);plot(log2(x),s(log2(x)),'color',[0 0 0]); hold on;
    
end
xlim([-9 1]);ylim([-7 0]);

%%
[x,y,s,a(2),b(2),a_d(2),a_u(2),b_d(2),b_u(2)] = SmoothHeatMap(A10,max(max(A10)),min(min(A10)),cutoff,gal,glu,fit_cuttoff,mid_value);
figure(2);plot(log2(x),log2(y),'o','color',[1, 165/255, 0],'markersize',1.5);hold on;
plot(log2(x),s(log2(x)),'color',[1, 165/255, 0])

[x,y,s,a(3),b(3),a_d(3),a_u(3),b_d(3),b_u(3)] = SmoothHeatMap(H3,max(max(H3)),min(min(H3)),cutoff,gal,glu,fit_cuttoff,mid_value);
figure(2);plot(log2(x),log2(y),'ob','markersize',1.5);hold on;
plot(log2(x),s(log2(x)),'b')

i=2
WT{i}(find(WT{i}==-inf))=nan;
WT{i}(find(WT{i}==0))=nan;

[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT{i},max(max(WT{i})),min(min(WT{i})),cutoff,gal_s288b_4283other,glu_s288b_4283other,fit_cuttoff,mid_value);
figure(2);
plot(log2(x),log2(y),'ok','markersize',1.5);
plot(log2(x),s(log2(x)),'k');hold on;

xlim([-9 2]);ylim([-9 0]);
legend({'s288c' 'A10' 'H3'})


%%
figure(3);
hold on;bar([1 2 3],a,'facecolor','none');
errorbar([1 2 3]',a',a'-a_d',a_u'-a','.k');ylim([0.65 1.05])
figure(4);hold on;bar([1 2 3],b,'facecolor','none');
errorbar([1 2 3]',b',b'-b_d',b_u'-b','.k');ylim([0.8 2.1])
Set_fig_YS(2,18,18,18)
Set_fig_YS(3,18,18,18)
Set_fig_YS(4,18,18,18)


%% differeent Dilution after 8 hrs.
close all
A10 = load('C:\Users\ys151\Dropbox\gal_paper\Data\20130806_dilution\output\heatmap_mean_A10');
H3 = load('C:\Users\ys151\Dropbox\gal_paper\Data\20130806_dilution\output\heatmap_mean_H3');


% Normalize with values of 1X8hrs values
A10_M{2} = flipud(A10.heatmaps1.Plate_1x_8H)

min_m = min(min(A10_M{2}));
max_m = max(max(A10_M{2}));
A10_M{1} = (flipud(A10.heatmaps1.Plate_1x)-min_m)/(max_m-min_m);
A10_M{2} = (flipud(A10.heatmaps1.Plate_1x_8H)-min_m)/(max_m-min_m);
A10_M{3} = (flipud(A10.heatmaps1.Plate_5_x)-min_m)/(max_m-min_m);
A10_M{4} = (flipud(A10.heatmaps1.Plate_5_x_8H)-min_m)/(max_m-min_m);
A10_M{5} = (flipud(A10.heatmaps1.Plate_fifth_x)-min_m)/(max_m-min_m);A10_M{5}(4,6) = nan;
A10_M{6} = (flipud(A10.heatmaps1.Plate_fifth_x_8H)-min_m)/(max_m-min_m);


H3_M{2} = flipud(H3.heatmaps1.Plate_1x_8H)

min_m = min(min(H3_M{2}));
max_m = max(max(H3_M{2}));
H3_M{1} = (flipud(H3.heatmaps1.Plate_1x)-min_m)/(max_m-min_m);
H3_M{2} = (flipud(H3.heatmaps1.Plate_1x_8H)-min_m)/(max_m-min_m);
H3_M{3} = (flipud(H3.heatmaps1.Plate_5_x)-min_m)/(max_m-min_m);
H3_M{4} = (flipud(H3.heatmaps1.Plate_5_x_8H)-min_m)/(max_m-min_m);
H3_M{5} = (flipud(H3.heatmaps1.Plate_fifth_x)-min_m)/(max_m-min_m);H3_M{5}(4,6) = nan;
H3_M{6} = (flipud(H3.heatmaps1.Plate_fifth_x_8H)-min_m)/(max_m-min_m);


gal = [0 2.^(-8:2)]
glu =[0 2.^(-6:0)];

for i = 1:6
    figure(1)
    subplot(3,2,i)
    h=pcolor(log2(gal),log2(glu),A10_M{i});hold on;contour(log2(gal),log2(glu),A10_M{i},cutoff,'k');axis square;
    set(h,'edgecolor','none')
    axis square;title('A10');
    xlim([-8 1]);ylim([-6 0]);
    colormap(cmap)
    
    figure(2)
    subplot(3,2,i)
    h=pcolor(log2(gal),log2(glu),H3_M{i});hold on;contour(log2(gal),log2(glu),H3_M{i},cutoff,'k');axis square;
    set(h,'edgecolor','none')
    axis square;title('H3');
    xlim([-8 1]);ylim([-6 0]);
    colormap(cmap)
    
end
Set_fig_YS(1,18,18,18);
Set_fig_YS(2,18,18,18);

%%
close all

figure(4);
title('A10')
hold on;
contour(log2(gal),log2(glu),A10_M{1},cutoff,'-k','linestyle','--');
contour(log2(gal),log2(glu),A10_M{2},cutoff,'k','linestyle','-','linewidth',2);
contour(log2(gal),log2(glu),A10_M{3},cutoff,'color',[1 165/255 0],'linestyle','--');
contour(log2(gal),log2(glu),A10_M{4},cutoff,'color',[1 165/255 0]);
contour(log2(gal),log2(glu),A10_M{5},cutoff,'b','linestyle','--');
contour(log2(gal),log2(glu),A10_M{6},cutoff,'b');

legend({'1X 4hrs' '1X 8hrs' '5X 4hrs' '5X 8hrs' '1/5X 4hrs' '1/5X 8hrs' });axis square

Set_fig_YS(4,18,18,18)

for i = 1:6
    
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(A10_M{i},1,0,cutoff,gal,glu,fit_cuttoff,mid_value);
    
end

figure(5)
ind = [5,6,1:4];
ind = [1 3 5 2 4 6]
% bar([1 2 4 5 7 8],b(ind),'facecolor','none');hold on;

ind = [1 2];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-ok');hold on
ind = [3 4];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-o','color',[1 165/255 0]);
ind = [5 6];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-ob');

title('A10')

Set_fig_YS(5,18,18,18);



figure(6);
title('H3')
hold on;
contour(log2(gal),log2(glu),H3_M{1},cutoff,'-k','linestyle','--');
contour(log2(gal),log2(glu),H3_M{2},cutoff,'k','linestyle','-','linewidth',2);
contour(log2(gal),log2(glu),H3_M{3},cutoff,'color',[1 165/255 0],'linestyle','--');
contour(log2(gal),log2(glu),H3_M{4},cutoff,'color',[1 165/255 0]);
contour(log2(gal),log2(glu),H3_M{5},cutoff,'b','linestyle','--');
contour(log2(gal),log2(glu),H3_M{6},cutoff,'b');

legend({'1X 4hrs' '1X 8hrs' '5X 4hrs' '5X 8hrs' '1/5X 4hrs' '1/5X 8hrs' });axis square

Set_fig_YS(6,18,18,18)

for i = 1:6
    
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(H3_M{i},1,0,cutoff,gal,glu,fit_cuttoff,mid_value);
    
end

figure(7)
% ind = [5,6,1:4];
% bar([1 2 4 5 7 8],b(ind),'facecolor','none');hold on;
% ind = [1 3 5 2 4 6]
% bar([1 2 4 5 7 8],b(ind),'facecolor','none');hold on;

ind = [1 2];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-ok');hold on
ind = [3 4];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-o','color',[1 165/255 0]);
ind = [5 6];errorbar([4 8]',b(ind)',b(ind)'-b_d(ind)',b_u(ind)'-b(ind)','-ob');

title('H3')
Set_fig_YS(7,18,18,18)


%%  Delitions
close all
clear wt
clear mt
clear a
clear b

gal = [0 2.^(-8:2)]
glu =[0 2.^(-6:0)];
mid_value = 4;
load('C:\Users\ys151\Dropbox\gal_paper\Data\Mutanta\Heatmap_matrices.mat');
name = {'Gal1' 'Gal2' 'Gal3' 'Gal4' 'Gal80'}

wt{1} = (heatmaps_wt.Gal1.mean.Matrix2);
wt{2} = (heatmaps_wt.Gal2.mean.Matrix2);
wt{3} = (heatmaps_wt.Gal3.mean.Matrix2);
wt{4} = (heatmaps_wt.Gal4.mean.Matrix2);
wt{5} = (heatmaps_wt.Gal80.mean.Matrix2);

mt{1} = (heatmaps_mut.Gal1.mean.Matrix1);
mt{2} = (heatmaps_mut.Gal2.mean.Matrix1);
mt{3} = (heatmaps_mut.Gal3.mean.Matrix1);
mt{4} = (heatmaps_mut.Gal4.mean.Matrix1);
mt{5} = (heatmaps_mut.Gal80.mean.Matrix1);

fit_cuttoff = [2^-1 2^-7];

for i = 1:5
    if i ==4
        fit_cuttoff = [2^-2.5 2^-7];
    end
    
    figure(1)
    wt{i}(find(wt{i}==-inf))=nan;
    min_wt = min(min(wt{i}));
    max_wt = max(max(wt{i}));
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(wt{i},max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);
    [x_mt,y_mt,s_mt,a_mt(i),b_mt(i),a_d_mt(i),a_u_mt(i),b_d_mt(i),b_u_mt(i)] = SmoothHeatMap(mt{i},max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);
    
    figure(1)
    subplot(5,1,i);
    plot(log2(x),log2(y),'.k');hold on;plot(log2(x),s(log2(x)),'k')
    plot(log2(x_mt),log2(y_mt),'.r');hold on;plot(log2(x_mt),s_mt(log2(x_mt)),'r')
    title(name{i});
    xlim([-8 1]);ylim([-6 0]);
    
    figure(2)
    subplot(5,1,i)
    bar([1 2],[a(i) a_mt(i)],'facecolor','none');hold on;
    errorbar([1 2 ]',[a(i) a_mt(i)]',[a(i) a_mt(i)]'-[a_d(i) a_d_mt(i)]',[a_u(i) a_u_mt(i)]'-[a(i) a_mt(i)]','.k')
    ylim([0.4 1.1])
    
    figure(3)
    subplot(5,1,i)
    bar([1 2],[b(i) b_mt(i)],'facecolor','none');hold on;
    errorbar([1 2 ]',[b(i) b_mt(i)]',[b(i) b_mt(i)]'-[b_d(i) b_d_mt(i)]',[b_u(i) b_u_mt(i)]'-[b(i) b_mt(i)]','.k')
    
    
end
Set_fig_YS(1,18,18,18)
Set_fig_YS(2,18,18,18)
Set_fig_YS(3,18,18,18)
%%
figure(4)
subplot(2,1,1)
bar([1:5],a./a_mt,'facecolor','none');refline(0,1);ylim([0.5 1.5])
subplot(2,1,2)
bar([1:5],b./b_mt,'facecolor','none');refline(0,1);ylim([0.5 1.5])
Set_fig_YS(4,18,18,18)

%% Compare GAl3 mig1 binding site delte with
clear b r
close all

cutoff = 0.2;

gal = [0 2.^(-8:2)]
glu =[0 2.^(-6:0)];
data = struct2cell(plates_hists.Second_plate_gal3_mut);
data = reshape(data,12,8);
for i = 1:12
    for j = 1:8
        wt_3_1(i,j) = data{i,j}.bfp_yfp.mean;
        mt_3_1(i,j) = data{i,j}.other_yfp.mean;
    end
end

max_wt = max(max(wt_3_1));
min_wt = min(min(mt_3_1));

i=1;
[x,y,s,a,b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(wt_3_1',max_wt,min_wt,cutoff,gal,glu);
[x_mt,y_mt,s_mt,a_mt(i),b_mt(i),a_d_mt(i),a_u_mt(i),b_d_mt(i),b_u_mt(i)] = SmoothHeatMap(mt_3_1',max_wt,min_wt,cutoff,gal,glu);



figure(1)
subplot(1,2,1);
h=pcolor(log2(gal),log2(glu),wt_3_1');
%     set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('WT diploid')
%     colormap(cmap)

subplot(1,2,2)
h=pcolor(log2(gal),log2(glu),mt_3_1');axis square;
%     set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('MT Haploid')

%     colormap(cmap)
figure(2)
plot(log2(x),log2(y),'k');hold on;plot(log2(x_mt),log2(y_mt),'r');
xlim([-8 1]);ylim([-6 0]);

%     min_wt = min(min(wt{3}));
%     max_wt = max(max(wt{3}));

[x,y,s,a,b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(wt{3},max_wt,min_wt,cutoff,gal,glu);
[x_mt,y_mt,s_mt,a_mt(i),b_mt(i),a_d_mt(i),a_u_mt(i),b_d_mt(i),b_u_mt(i)] = SmoothHeatMap(mt{3},max_wt,min_wt,cutoff,gal,glu);

plot(log2(x),log2(y),'--k');hold on;plot(log2(x_mt),log2(y_mt),'--r');

legend({'WT-D','MT: Gal3-Mig1' 'WT-D' 'MT: hetro\Delta'})
Set_fig_YS(2,18,18,18)

%% large gradiant GAl3-Mig1 delete

load('C:\Users\ys151\Dropbox\Gal3mig1delete\LargeGradiaent\plates_hists_EMD')
clear b r wt mt x
close all
cmap =  buildcmap_YS([0 0 0.5;1 1 1;0.5 0 0]);

cutoff = 0.5;

gal = [0 2.^(-8:2)]
glu =[0 2.^(-6:0)];

data = struct2cell(plates_hists);

for k = 1:6
    if k==6
        xx = reshape( struct2cell(data{k}),12,7)
        for i = 1:12
            for j = 1:7
                wt(i,j,k) = xx{i,j}.bfp_yfp.mean;
                mt(i,j,k) = xx{i,j}.other_yfp.mean;
            end
        end
    else
        xx = reshape( struct2cell(data{k}),12,8)
        
        for i = 1:12
            for j = 1:8
                wt(i,j,k) = xx{i,j}.bfp_yfp.mean;
                mt(i,j,k) = xx{i,j}.other_yfp.mean;
            end
        end
    end
end

W = [ wt(:,:,1)' wt(:,:,2)' ; wt(:,:,3)' wt(:,:,4)';nan*wt(:,:,5)' wt(:,:,6)'];
W(W==0)=nan;
W(W==-inf)=nan
Mt = [ mt(:,:,1)' mt(:,:,2)' ; mt(:,:,3)' mt(:,:,4)';mt(:,:,5)' mt(:,:,6)'];
Mt(Mt==0)=nan;
Mt(Mt==-inf)=nan;

gal = [0 2.^[-9:0.5:2]]
glu = [0 0 2.^[-9.5:0.5:1]]

fit_cuttoff = [2^-1 2^-7];
mid_value = 2^-4;

fit_cuttoff = 2.^[-1 -8];

figure(1)
subplot(1,2,1)
h=pcolor(log2(gal),log2(glu),W);hold on;
axis square;
set(h,'edgecolor','none')
axis square;title('WT');
xlim([-9 1]);ylim([-9 0]);
colormap(cmap)
[x,y,s,a(1),b(1),a_d(1),a_u(1),b_d(1),b_u(1)] = SmoothHeatMap(W,max(max(W)),min(min(W)),cutoff,gal,glu,fit_cuttoff,mid_value);
plot(log2(x),log2(y),'k');
figure(2)
plot(log2(x),log2(y),'.k');hold on;plot(log2(x),s(log2(x)),'k');hold on

figure(1);subplot(1,2,2)
h=pcolor(log2(gal),log2(glu),Mt);hold on;
axis square;
set(h,'edgecolor','none')
axis square;title('MT');
xlim([-9 1]);ylim([-9 0]);
colormap(cmap)
[x,y,s,a(2),b(2),a_d(2),a_u(2),b_d(2),b_u(2)] = SmoothHeatMap(Mt,max(max(W)),min(min(W)),cutoff,gal,glu,fit_cuttoff,mid_value);
plot(log2(x),log2(y),'k');
figure(2)
plot(log2(x),log2(y),'.r');hold on;plot(log2(x),s(log2(x)),'r');hold on;title('Black = WT; Red = MT');
xlim([-9 1]);ylim([-9 0]);



figure(3)
bar([1 2],[a],'facecolor','none');hold on;
errorbar([1 2 ]',a,a-a_d,a_u-a,'.k')
xlim([0 3]);ylim([0.5 1.5]);title('Slope: left:WT right: MT')

figure(4)
bar([1 2],[b],'facecolor','none');hold on;
errorbar([1 2 ]',b,b-b_d,b_u-b,'.k')
xlim([0 3]);ylim([0.5 2]);title('Cutoff: left:WT right: MT')

Set_fig_YS(1,18,18,18)
Set_fig_YS(2,18,18,18)
Set_fig_YS(3,18,18,18)
Set_fig_YS(4,18,18,18)
%% Gal4 Mig1 delte
load('C:\Users\ys151\Dropbox\Gal4Mig1Delete\plates_hists_EMD')

clear b r
close all

cmap =  buildcmap_YS([1 1 1;0 0 0.5;0 0.5 0;0.5 0 0]);


cutoff = 0.2;
fit_cuttoff = 2.^[-1 -7];
mid_value = 2^-4;

gal = [0 2.^(-8:2)]
glu =[0 2.^(-6:0)];

data = struct2cell(plates_hists.YM4277);
data = reshape(data,12,8);
for i = 1:12
    for j = 1:8
        wt_s288c_77(i,j) = data{i,j}.bfp_yfp.mean;
        wt_77(i,j) = data{i,j}.other_yfp.mean;
    end
end

data = struct2cell(plates_hists.YM4283);
data = reshape(data,12,8);
for i = 1:12
    for j = 1:8
        wt_s288c_83(i,j) = data{i,j}.bfp_yfp.mean;
        wt_83(i,j) = data{i,j}.other_yfp.mean;
    end
end

wt_s288c_83(wt_s288c_83==-inf)=nan;
max_wt = max(max(wt_s288c_83));
min_wt = min(min(wt_s288c_83));

[x1,y1,s,a(1),b(1),a_d(1),a_u(1),b_d(1),b_u(1)] = SmoothHeatMap(wt_s288c_83',max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);
[x2,y2,s,a(2),b(2),a_d(2),a_u(2),b_d(2),b_u(2)] = SmoothHeatMap(wt_83',max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);
[x3,y3,s,a(3),b(3),a_d(3),a_u(3),b_d(3),b_u(3)] = SmoothHeatMap(wt_s288c_77',max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);
[x4,y4,s,a(4),b(4),a_d(4),a_u(4),b_d(4),b_u(4)] = SmoothHeatMap(wt_77',max_wt,min_wt,cutoff,gal,glu,fit_cuttoff,mid_value);


wt_s288c_83 = mat2gray(wt_s288c_83, [min_wt max_wt])
wt_83 = mat2gray(wt_83, [min_wt max_wt])
wt_s288c_77 = mat2gray(wt_s288c_77, [min_wt max_wt])
wt_77 = mat2gray(wt_77, [min_wt max_wt])

figure(1)
subplot(2,2,1);
h=pcolor(log2(gal),log2(glu),77');
set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('77');
colormap(cmap);

subplot(2,2,2);
h=pcolor(log2(gal),log2(glu),wt_83');
set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('83');
colormap(cmap);

subplot(2,2,3);
h=pcolor(log2(gal),log2(glu),wt_s288c_77');
set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('s288c (with 77)');
colormap(cmap);

subplot(2,2,4);
h=pcolor(log2(gal),log2(glu),wt_s288c_83');
set(h,'edgecolor','none')
axis square;
xlim([-8 1]);ylim([-6 0]);
title('s288c (with 83)');
colormap(cmap);

Set_fig_YS(1,18,18,18)


figure(2)
plot(log2(x1),log2(y1),'k');hold on;plot(log2(x2),log2(y2),'r');h
xlim([-8 1]);ylim([-6 0]);

legend({'s288c','4283'})
Set_fig_YS(2,18,18,18)