function [M,D]=plot_S_plots_1C_S6()

[D,D2] = load_data_S_plots()
save('Data_S_plots','D','D2')

%Double gradient May 17 2014
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

filename='output/Figure_1C_S_plot_WT_20150729';
export_fig(filename, '-pdf','-transparent','-nocrop');

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

end
