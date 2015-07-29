function plot_figure_S6C_S_plot

%Projection of the response to the galactose/glucose axis, as in Fig. 1C. 
%Note that the mean projection does not collapse as tightly onto a single curve as the other two metrics

load('../data/20140701_stitched_areas/output/Data_S_plots.mat');

%% Make the S-plots for different mesures, remvove points -7
color_vec = [0 0 0;1 0 0;0 0 1;0 1/2 0];

%Double gradient May 17 2014
glc_final = [0 2.^[-7:0.5:0]];
gal_final = [0 2.^[-9:0.5:2]];


%close all
gal= gal_final(2:end);
glc = glc_final(2:end);

clear rat
for i = 1:length(glc)
    for j = 1:length(gal)
        
        rat(i,j) = gal(j)/glc(i);
        
    end
end


figure(7)
for i = 1:3
    %     figure(i)
    
    plot(log10(rat),D2{i},'o','color',color_vec(i,:),'markersize',4,'markerfacecolor',color_vec(i,:));hold on;
    
    [XOut, YOut, ZOut] = prepareSurfaceData(gal, glc, D2{i});
    s = fit([XOut, YOut], ZOut,'c./(c + (y./x).^n)');
end

box off;
ylim([0 1]);

ylabel('Fraction of cells inducing')
xlabel('log_{10} (Galactose/Glucose) ')
Set_fig_RE(figure(7),18,18,22);


filename='Figure_S6C_S_plot_WT_different_metrics';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');
