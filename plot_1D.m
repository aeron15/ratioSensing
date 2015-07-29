function plot_1D()
%PLOT_1D
%Comparison of models of signal integration (SI Appendix) by threshold sensing (Upper) and ratio sensing (Lower), displayed as in B and C.
%Inspired by S288cSurfaceFit

%% Pure ratio.  Heatmap and S plot.
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
colormap(cmap);caxis = [0 1];
box off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);

for i = 1:length(b)
    for j = 1:length(a)
        
        rat(i,j) = a(j)/b(i);
        
    end
end
subplot(1,2,2)

plot(rat,ff,'color','k','linewidth',2);hold on;
set(gca,'xscale','log');

box off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([10^-4,10^4]);ylim([-0.01 1.01])
axis('square');

Set_fig_RE(figure(10),12,18,18)

filename='Figure1D_heatmap_S_plot_ratio_sensor';
export_fig_specific_path(filename, '-pdf','-nocrop');

%% Threshold model. Heatmap and S plot
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
colormap(cmap);caxis = [0 1];
box off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);

for i = 1:length(b)
    for j = 1:length(a)
        
        rat(i,j) = a(j)/b(i);
        
    end
end
subplot(1,2,2)

for i = 1:length(c)
    x =  logspace(-4,5,100);
    ff = (1./(1 + (1./(sqrt(c(i))*sqrt(x))).^2)).*(1./(1 + ((sqrt(c(i))./sqrt(x))./1).^2));
    plot(x,ff,'color','k','linewidth',2);hold on;
end

set(gca,'xscale','log');

box off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([10^-4,10^5]);ylim([-0.01 1.01])
axis('square');

Set_fig_RE(figure(10),12,18,18)

filename='Figure1D_heatmap_S_plot_threshold_sensor';
export_fig_specific_path(filename, '-pdf','-nocrop');

end

