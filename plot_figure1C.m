function plot_figure1C()
% Plot figure 1C.Fraction of inducing cells as a function of the ratio of galactose and glucose concentrations.
load('../data/20140701_stitched_areas/output/Data_S_plots.mat');

glc_final = [0 2.^[-9:0.5:0]];
gal_final = [0 2.^[-9:0.5:2]];

gal= gal_final(2:end);
glc = glc_final(2:end);

%%
for i = 1:length(glc)
    for j = 1:length(gal)
        
        rat(i,j) = gal(j)/glc(i);
        
    end
end

color_vec = [0 0 0];

figure(6)
for i = 1
    plot(rat,D{i},'o','color',color_vec(i,:),'markersize',4,'markerfacecolor',color_vec(i,:));hold on;
    
    set(gca,'xscale','log');
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc,D{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
    plot(rat,s.c./(s.c + (rat.^-s.n)),'color',color_vec(i,:),'linewidth',4)
end

box off;xlim([10^-3,10^3.4]);ylim([0 1]);

ylabel('Fraction of cells inducing')
xlabel('log_{10} (Galactose/Glucose) ')

Set_fig_RE(figure(6),14,14,14);

filename='Figure_1C_S_plot_WT';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');



