function plot_Figure4B()
%PLOT_FIGURE4B
%ratio output can be generally achieved by a simple module in which two input molecules, 
%an activator (green) and a repressor (blue), bind to an integrator molecule?a promoter, 
%transporter, scaffold protein, etc. Mutual inhibition, e, is necessary for a robust ratio response. 
%% Figure. Models of ratio sensing.

cmap =  buildcmap_YS([1 1 1;0.5 0 0]);
% cmap =  buildcmap('br');
cmap=cbrewer('seq', 'YlOrRd', 15);

a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);


aa = logspace(-2,10,1000000);
bb = logspace(-2,10,1000000);

d = [10^0 10^3 10^7];
Kta = [10^10 10^2 10^2];
Ktb = [10^10 10^10 10^2];

figure(1)

for i = 1:length(d)
    for k = 1%length(Kta)
        a = logspace(-2,5,50);
        b = logspace(-2,5,50);
        [A,B] = meshgrid(a,b);
        
        B = B*Ktb(k)./(Ktb(k) + B);
        A = A*Kta(k)./(Kta(k) + A);
        ff = (1./(1 + B./A +1./A + B/d(i)));
        
        subplot(1,3,3*(k-1)+i)
        h = contourf(a,b,ff,5);hold on;
        
        %To highlight the
        h = surf(a,b,ff);hold on;
        
        plot3(aa,1./(1-aa/d(i)),ones(size(aa))*10^5,'color',[0 0 0],'linewidth',2);
        plot3(ones(size(a))*Kta(k),b,ones(size(a))*10^5,'--k','linewidth',1);
        plot3(a,ones(size(a))*Ktb(k),ones(size(a))*10^5,'--k','linewidth',1);
        set(h,'edgecolor','none','facecolor','none');
        set(gca,'xscale','log','yscale','log','zscale','linear');
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
        %         xlabel('A/K_{A}');ylabel('R/K_{R}');
        colormap(cmap);
        caxis([0 1]);
        colorbar('location', 'NorthOutside')
    end
end


Set_fig_RE(figure(1),12,12,16)
filename=['Figure4B_Models_signal_integration'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

end

