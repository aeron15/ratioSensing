%% Analyze different metrics s288c histograms.


clear all
close all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load histograms and concentrations

% s288c replicae 1
%load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\\20131019_2D_LG\output\plates_hists_EMD')
load ('../data/s288c/20131019_2D_LG/output/plates_hists_EMD')
data = struct2cell(plates_hists);
th_const = 2.5;
off_peak = 2;
for k = 1:4
    clear mean_temp prec_temp prec_area_temp
    

    xx = struct2cell(data{k});
    for i = 1:length(xx)
        if k == 1 && i==1
            y_off = xx{1}.bfp_yfp.f;
            x_off = xx{1}.bfp_yfp.xi;
        end
        [prec_const,prec_adaptive,prec_area_temp(i)] = AnalyzeHist(xx{i}.bfp_yfp.f,xx{i}.bfp_yfp.xi,th_const,off_peak,y_off,x_off);
        mean_temp(i) = xx{i}.bfp_yfp.mean;
        prec_temp(i) = xx{i}.bfp_yfp.perc_ind;
    end
    if k==3
        mean_temp(57:end+1) = mean_temp(56:end);mean_temp(56)=nan;
        prec_temp(57:end+1) = prec_temp(56:end);prec_temp(56)=nan;
        prec_area_temp(57:end+1) = prec_area_temp(56:end);prec_area_temp(56)=nan;

    end
    mean_mat(:,:,k) =  (reshape(mean_temp,12,8));
    prec_mat(:,:,k) =  (reshape(prec_temp,12,8));
    prec_area_mat(:,:,k) =  (reshape(prec_area_temp,12,8));

end
    
WT_mean{1} = [ mean_mat(:,:,1)',mean_mat(:,:,2)' ; mean_mat(:,:,3)' ,mean_mat(:,:,4)'];
WT_prec{1} = [ prec_mat(:,:,1)',prec_mat(:,:,2)' ; prec_mat(:,:,3)' ,prec_mat(:,:,4)'];
WT_prec_area{1} = [ prec_area_mat(:,:,1)',prec_area_mat(:,:,2)' ; prec_area_mat(:,:,3)' ,prec_area_mat(:,:,4)'];

%% s288c replicae 2

load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\20140217_4283\output\plates_hists_EMD')
data = struct2cell(plates_hists);

for k = 1:6
    clear mean_temp prec_temp prec_area_temp
    
    xx = struct2cell(data{k});
    for i = 1:length(xx)
        
        if k == 1 && i==1
            xxx = struct2cell(data{5});

            y_off = xx{1}.bfp_yfp.f;
            x_off = xx{1}.bfp_yfp.xi;
        end
        [prec_const,prec_adaptive,prec_area_temp(i)] = AnalyzeHist(xx{i}.bfp_yfp.f,xx{i}.bfp_yfp.xi,th_const,off_peak,y_off,x_off);
        
        mean_temp(i) = xx{i}.bfp_yfp.mean;
        prec_temp(i) = xx{i}.bfp_yfp.perc_ind;

    end
    if k<=4
        mean_mat(:,:,k) =  (reshape(mean_temp,12,8));
        prec_mat(:,:,k) =  (reshape(prec_temp,12,8));
        prec_area_mat(:,:,k) =  (reshape(prec_area_temp,12,8));

    else
        mean_mat(:,:,k) =  [reshape(mean_temp,12,4),nan*ones(12,4)];
        prec_mat(:,:,k) =  [reshape(prec_temp,12,4),nan*ones(12,4)];
        prec_area_mat(:,:,k) =[reshape(prec_area_temp,12,4),nan*ones(12,4)];
    end
    
    
end
WT_mean{2} = [ mean_mat(:,:,1)',mean_mat(:,:,2)' ; mean_mat(:,:,3)' ,mean_mat(:,:,4)';mean_mat(:,1:4,5)' ,mean_mat(:,1:4,6)'];
WT_prec{2} = [ prec_mat(:,:,1)',prec_mat(:,:,2)' ; prec_mat(:,:,3)' ,prec_mat(:,:,4)';prec_mat(:,1:4,5)' ,prec_mat(:,1:4,6)'];
WT_prec_area{2} = [ prec_area_mat(:,:,1)',prec_area_mat(:,:,2)' ; prec_area_mat(:,:,3)' ,prec_area_mat(:,:,4)';prec_area_mat(:,1:4,5)' ,prec_area_mat(:,1:4,6)'];


%% s288c replicae 3

load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\\20140217_YM4277_mig1_del_Gal4\output\plates_hists_EMD')
data = struct2cell(plates_hists);

for k = 1:6
    clear mean_temp prec_temp prec_area_temp
    
    xx = struct2cell(data{k});
    for i = 1:length(xx)
        
        if k == 1 && i==1
            xxx = struct2cell(data{5});

            y_off = xx{1}.bfp_yfp.f;
            x_off = xx{1}.bfp_yfp.xi;
        end
        [prec_const,prec_adaptive,prec_area_temp(i)] = AnalyzeHist(xx{i}.bfp_yfp.f,xx{i}.bfp_yfp.xi,th_const,off_peak,y_off,x_off);
        
        mean_temp(i) = xx{i}.bfp_yfp.mean;
        prec_temp(i) = xx{i}.bfp_yfp.perc_ind;

    end
    if k<=4
        if k==2
            mean_temp(3:end+1) = mean_temp(2:end);mean_temp(2)=nan;
            prec_temp(3:end+1) = prec_temp(2:end);prec_temp(2)=nan;
            prec_area_temp(3:end+1) = prec_area_temp(2:end);prec_area_temp(2)=nan;
            
        end
    
        mean_mat(:,:,k) =  (reshape(mean_temp,12,8));
        prec_mat(:,:,k) =  (reshape(prec_temp,12,8));
        prec_area_mat(:,:,k) =  (reshape(prec_area_temp,12,8));

    else
        mean_mat(:,:,k) =  [reshape(mean_temp,12,4),nan*ones(12,4)];
        prec_mat(:,:,k) =  [reshape(prec_temp,12,4),nan*ones(12,4)];
        prec_area_mat(:,:,k) =[reshape(prec_area_temp,12,4),nan*ones(12,4)];
    end

    
end
WT_mean{3} = [ mean_mat(:,:,1)',mean_mat(:,:,2)' ; mean_mat(:,:,3)' ,mean_mat(:,:,4)';mean_mat(:,1:4,6)' ,mean_mat(:,1:4,5)'];
WT_prec{3} = [ prec_mat(:,:,1)',prec_mat(:,:,2)' ; prec_mat(:,:,3)' ,prec_mat(:,:,4)';prec_mat(:,1:4,6)' ,prec_mat(:,1:4,5)'];
WT_prec_area{3} = [ prec_area_mat(:,:,1)',prec_area_mat(:,:,2)' ; prec_area_mat(:,:,3)' ,prec_area_mat(:,:,4)';prec_area_mat(:,1:4,6)' ,prec_area_mat(:,1:4,5)'];


gal{1} = [0 2.^[-9:0.5:2]];
glu{1} = [0 2.^[-7:0.5:0]];

gal{2} = [0 2.^[-9:0.5:2]];
glu{2} = [0 2.^[-9:0.5:0]];

gal{3} = [0 2.^[-9:0.5:2]];
glu{3} = [0 2.^[-9:0.5:0]];


%% PLOT the MAIN replicates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all


n=3
figure(1)
for i = 1:n
    subplot(1,3,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_mean{i});
    axis square;xlim([-9 1]);ylim([-7 0]);title('Normalized Mean')

end

figure(2)
for i = 1:n
    subplot(1,3,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_prec{i});
    axis square;xlim([-9 1]);ylim([-7 0]);title('Constant threshold')

end

figure(3)
for i = 1:n
    subplot(1,3,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_prec_area{i});
    axis square;xlim([-9 1]);ylim([-7 0]);title('area')

end

% plot the Mean decision fronts

fit_cuttoff = [2^-2 2^-6];
mid_value = 2^-4;
cutoff=0.2;



for i = 1:n
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = min(min(WT_mean{i}));
    max_wt = max(max(WT_mean{i}));
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_mean{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(4);subplot(2,2,[1 3]);plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 0],'markersize',1.5);hold on
    plot(log2(x),s(log2(x)),'color',[0 0 0]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]);

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.k');ylim([0 2])
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.k');ylim([0.4 2.1])

figure(4);title('Normalized Mean');Set_fig_YS(figure(4),18,18,18);

% plot the prec decision fronts



for i = 1:n
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(5);subplot(2,2,[1 3]);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[1 0 0],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[1 0 0]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]); title('Constant Treshold ');

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.r');ylim([0 2]);
% title(['Slope = ',num2str(
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.r');ylim([0.4 2.1])

Set_fig_YS(figure(5),18,18,18);


% plot the prec area decision fronts



for i = 1:n
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec_area{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(6);subplot(2,2,[1 3]);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 1],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[0 0 1]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]);

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.r');ylim([0 2])
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.r');ylim([0.4 2.1])
figure(6);title('Area');Set_fig_YS(figure(6),18,18,18);

% Compare mean to prec

for i = 1:3
    
    
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = min(min(WT_mean{i}));
    max_wt = max(max(WT_mean{i}));
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_mean{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(1,3,i);plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 0],'markersize',1.5);hold on
    plot(log2(x),s(log2(x)),'color',[0 0 0]); hold on;
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(1,3,i);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[1 0 0],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[1 0 0]); hold on;
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec_area{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(1,3,i);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 1],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[0 0 1]); hold on;
    
    
    
    
    
    axis square;xlim([-9 1]);ylim([-7 0]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Include 96 Well plates replucates.
close all
load('C:\Users\ys151\Dropbox\gal_paper\Data\Delitions\20131114_galdeletions\output\plates_hists_EMD');
data = struct2cell(plates_hists);

th_const = 2.5;
off_peak = 2;

for k = 1:5
     clear mean_temp prec_temp prec_area_temp
    
        xx = struct2cell(data{k});
        for i = 1:length(xx)
            
            if i==1
                xxx = struct2cell(data{k});
                
                y_off = xx{1}.mCh_yfp.f;
                x_off = xx{1}.mCh_yfp.xi;
            end
            [prec_const,prec_adaptive,prec_area_temp(i)] = AnalyzeHist(xx{i}.mCh_yfp.f,xx{i}.mCh_yfp.xi,th_const,off_peak,y_off,x_off);
            
            mean_temp(i) = xx{i}.mCh_yfp.mean;
            prec_temp(i) = xx{i}.mCh_yfp.perc_ind;
            
        end
        
        mean_mat(:,:,k) =  (reshape(mean_temp,12,8));
        prec_mat(:,:,k) =  (reshape(prec_temp,12,8));
        prec_area_mat(:,:,k) =  (reshape(prec_area_temp,12,8));
        
        
    WT_mean{k+3}  = mean_mat(:,:,k)';
    WT_prec{k+3} = prec_mat(:,:,k)' ;
    WT_prec_area{k+3} =prec_area_mat(:,:,k)' ;
end
%%
load('C:\Users\ys151\Dropbox\gal_paper\Data\Delitions\20140227_DG_gal_deletions\output\plates_hists_EMD');
data = struct2cell(plates_hists);

th_const = 2.5;
off_peak = 2;

for k = 1:5
     clear mean_temp prec_temp prec_area_temp
    
        xx = struct2cell(data{k});
        for i = 1:length(xx)
            
            if i==1
                xxx = struct2cell(data{k});
                
                y_off = xx{1}.bfp_yfp.f;
                x_off = xx{1}.bfp_yfp.xi;
            end
%             [prec_const,prec_adaptive,prec_area_temp(i)] = AnalyzeHist(xx{i}.bfp_yfp.f,xx{i}.bfp_yfp.xi,th_const,off_peak,y_off,x_off);
            
            mean_temp(i) = xx{i}.bfp_yfp.mean;
            prec_temp(i) = xx{i}.bfp_yfp.perc_ind;
            
        end
        
        mean_mat(:,:,k) =  (reshape(mean_temp,12,8));
        prec_mat(:,:,k) =  (reshape(prec_temp,12,8));
%         prec_area_mat(:,:,k) =  (reshape(prec_area_temp,12,8));
        
        
    WT_mean{k+8}  = mean_mat(:,:,k)';
    WT_prec{k+8} = prec_mat(:,:,k)' ;
    WT_prec_area{k+8} =prec_area_mat(:,:,k)' ;
end

% WT_mean{3} = (heatmaps_wt.Gal1.mean.Matrix2);
% WT_mean{4} = (heatmaps_wt.Gal2.mean.Matrix2);
% WT_mean{5} = (heatmaps_wt.Gal3.mean.Matrix2);
% WT_mean{6} = (heatmaps_wt.Gal4.mean.Matrix2);
% WT_mean{7} = (heatmaps_wt.Gal80.mean.Matrix2);
% 
% WT_prec{3} = (heatmaps_wt.Gal1.perc_ind.Matrix2);
% WT_prec{4} = (heatmaps_wt.Gal2.perc_ind.Matrix2);
% WT_prec{5} = (heatmaps_wt.Gal3.perc_ind.Matrix2);
% WT_prec{6} = (heatmaps_wt.Gal4.perc_ind.Matrix2);
% WT_prec{7} = (heatmaps_wt.Gal80.perc_ind.Matrix2);

%%
gal{1} = [0 2.^[-9:0.5:2]];
glu{1} = [0 2.^[-7:0.5:0]];

gal{2} = [0 2.^[-9:0.5:2]];
glu{2} = [0 2.^[-9:0.5:0]];

for i = 4:9
    gal{i} = [0 2.^(-8:2)];
    glu{i} =[0 2.^(-7:-1)];
    
end
for i = 9:13
    gal{i} = [0 2.^(-8:2)];
    glu{i} =[0 2.^(-6:0)];
    
end
%% PLOT the mean heat map

n = 8;

figure(1)
for i = 1:n
    subplot(4,4,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_mean{i});
    axis square;xlim([-9 1]);ylim([-7 0]);

end

figure(2)
for i = 1:n
    subplot(4,4,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_prec{i});
    axis square;xlim([-9 1]);ylim([-7 0]);

end

figure(3)
for i = 1:n
    subplot(4,4,i)
    pcolor(log2(gal{i}),log2(glu{i}),WT_prec_area{i});
    axis square;xlim([-9 1]);ylim([-7 0]);

end


%% plot the Mean decision fronts

fit_cuttoff = [2^-2 2^-6];
mid_value = 2^-4;
cutoff=0.2;



for i = 1:n
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = min(min(WT_mean{i}));
    max_wt = max(max(WT_mean{i}));
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_mean{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(4);subplot(2,2,[1 3]);plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 0],'markersize',1.5);hold on
    plot(log2(x),s(log2(x)),'color',[0 0 0]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]);

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.k');ylim([0 2])
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.k');ylim([0.4 2.1])


%% plot the prec decision fronts



for i = 1:n
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(5);subplot(2,2,[1 3]);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[1 0 0],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[1 0 0]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]);

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.r');ylim([0 2])
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.r');ylim([0.4 2.1])


%% plot the prec area decision fronts



for i = 1:n
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec_area{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(6);subplot(2,2,[1 3]);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 1],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[0 0 1]); hold on;
    
end
axis square;xlim([-9 1]);ylim([-7 0]);

subplot(2,2,2)
hold on;bar([1:n],a,'facecolor','none');
errorbar([1:n]',a',a'-a_d',a_u'-a','.r');ylim([0 2])
subplot(2,2,4)
hold on;bar([1:n],b,'facecolor','none');
errorbar([1:n]',b',b'-b_d',b_u'-b','.r');ylim([0.4 2.1])

%% Compare mean to prec

for i = 1:n
    
    
    
    WT_mean{i}(find(WT_mean{i}==-inf))=nan;
    min_wt = min(min(WT_mean{i}));
    max_wt = max(max(WT_mean{i}));
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_mean{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(3,3,i);plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 0],'markersize',1.5);hold on
    plot(log2(x),s(log2(x)),'color',[0 0 0]); hold on;
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(3,3,i);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[1 0 0],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[1 0 0]); hold on;
    
    min_wt = 0;
    max_wt = 1;
    [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT_prec_area{i},max_wt,min_wt,cutoff,gal{i},glu{i},fit_cuttoff,mid_value);
    figure(7);subplot(3,3,i);
    plot(log2(x),log2(y),'o','markerfacecolor',[1 1 1],'color',[0 0 1],'markersize',1.5);hold on;
    plot(log2(x),s(log2(x)),'color',[0 0 1]); hold on;
    
    
    
    
    axis square;xlim([-9 1]);ylim([-7 0]);

end
