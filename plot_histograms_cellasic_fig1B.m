function plot_histograms_cellasic_fig1B()
%PLOT_HISTOGRAMS_CELLASIC_FIG1B
%Comparison of cells monitored by live mi- croscopy to FCM at three sugar mixtures, denoted by numbered squares
% inspired by CellAcisMain_17_APR_2014_Gal1_Dynamics which compiles all the
% data from the server

close all;

load('../data/Gal1DynamicsOnix_17_APR_2014/data.mat')

dt=15/60;
T =[0:49]*dt;

%%
ctrs = linspace(80,1000,100);
for i = 1:4
    for j=1:50
        y = yfp_hist{5*(i-1)+1:5*(i-1)+5,j};
        yfp_hist_conc(i,j,:) = hist(y,ctrs);
        yfp_hist_conc(i,j,:) = yfp_hist_conc(i,j,:) /max(yfp_hist_conc(i,j,:) );
        med(i,j) = median(y);
        
        y = rfp_hist{5*(i-1)+1:5*(i-1)+5,j};
        rfp_hist_conc(i,j,:) = hist(y,ctrs);
        rfp_hist_conc(i,j,:) = rfp_hist_conc(i,j,:) /max(rfp_hist_conc(i,j,:) );
        med_r(i,j) = median(y);
        
    end
end

%% Plot histograms microscopy

time_probe = [32]+1;

TT = T(time_probe);

k =1;
for i = 1:3
    x = squeeze(yfp_hist_conc(i,time_probe,:));
    
    x=yfp_hist_conc(i,time_probe,:);
    
    figure(3)
    for j = 1:length(TT)
        subplot(4,length(time_probe),k)
        plot(ctrs,x(j,:),'k','linewidth',2);
        set(gca,'xscale','log');xlim([80,1000]);ylim([0 1.1]);
        axis square
        box off
        set(gcf,'color','None')
        set(gca,'YTick',[],'YTickLabel','')
        k = k+1;
        
    end
end

filename='Figure_1B_microscopy_pannels';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

end