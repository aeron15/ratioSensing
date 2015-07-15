function figure2

%% Figure. Models of ratio sensing.
close all
%clear all

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
        a = logspace(-2,5,100);
        b = logspace(-2,5,100);
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
%         set(h,'edgecolor','none','facecolor','interp');
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
filename=['Models_signal_integration'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
%% transport
d = [10^0 10^3 10^7];
Kta = [10^10 10^2 10^2];
Ktb = [10^10 10^10 10^2];


for i = 1:length(d)
    for k = 1:length(Kta)
        a = logspace(-2,5,100);
        b = logspace(-2,5,100);
        [A,B] = meshgrid(a,b);
        
        B = B*Ktb(k)./(Ktb(k) + B);
        A = A*Kta(k)./(Kta(k) + A);
        ff = (1./(1 + B./A +1./A + B/d(i)));
        
        figure(2);
        subplot(3,3,3*(k-1)+i);
        h = contourf(a,b,ff,5);hold on;
        plot3(aa,1./(1-aa/d(i)),ones(size(aa))*10^5,'color',[0 0 0],'linewidth',1);
        plot3(ones(size(a))*Kta(k),b,ones(size(a))*10^5,'--k','linewidth',1);
        plot3(a,ones(size(a))*Ktb(k),ones(size(a))*10^5,'--k','linewidth',1);
%         set(h,'edgecolor','none','facecolor','interp');
        set(gca,'xscale','log','yscale','log','zscale','linear');
        set(gca,'xtick',[1  100   10000]);        set(gca,'xticklabel',[1  2   4]);
        set(gca,'ytick',[1  100   10000]);set(gca,'yticklabel',[1  2   4]);
        
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
%         xlabel('A/K_{A}');ylabel('R/K_{R}');
        colormap(cmap);caxis([0 1]);
    
    end
end
Set_fig_YS(figure(2),9,22,22)


%% activator is a ratio 
close all
clear all
cmap=cbrewer('seq', 'YlOrRd', 25);

% cmap =  buildcmap_YS([1 1 1;0.5 0 0]);
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);


aa = logspace(-2,5,10000);
bb = logspace(-2,5,10000);

d = [1 10^2 10^20];
% Kta = [10^10 10^2 10^2];
% Ktb = [10^10 10^10 10^2];

figure(2)
KA = 10;
KR = 100;

CA = 0.01;
rho = 1;
for i = 1:length(d)
    for k = 1
        a = logspace(-2,5,100);
        b = logspace(-2,5,100);
        [A,B] = meshgrid(a,b);
                aa = logspace(-2,5,10000);

        %     B = B*Ktb(k)./(Ktb(k) + B);
        %     A = 1./(1+B./A).*10^5;
        %
        
        %     B = 1000*B./(B + 1000);
        ff = 1./(1 + CA*rho + B*(1/d(i) + CA*rho) + (B./A).*(B+1));
        
        subplot(1,3,3*(k-1)+i)
        h = contourf(a,b,ff,5);hold on;
        %        plot3(rho*b,b,ones(size(a))*10^5,'k','linewidth',1);
        %
%         plot3(a,a*rho*CA,ones(size(a))*10^5,'k','linewidth',1);
        
        plot3( aa,aa*(1/d(i)+rho*CA)-1,ones(size(aa))*10^5,'k','linewidth',1);
        plot(a,ones(size(a)),'k','linewidth',1)
%         set(h,'edgecolor','none','facecolor','interp');
        set(gca,'xscale','log','yscale','log','zscale','linear');
            set(gca,'xtick',[1  100   10000]);        set(gca,'xticklabel',[1  2   4]);
        set(gca,'ytick',[1  100   10000]);set(gca,'yticklabel',[1  2   4]);
        axis('square');view([0 90]);
        xlim([10^-2 10^5]);ylim([10^-2 10^5])
%         xlabel('s_A/C_{A}');ylabel('s_R/C_{R}');
        colormap(cmap);caxis([0 1]);
    end
end
box off
grid off
Set_fig_YS(figure(2),9,9,9)
% figure(3)
% subplot(2,1,1)
% h = surf(a,b,A);hold on;
% set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','linear');
% axis('square');view([0 90]);
% xlim([0 10^5]);ylim([0 10^5])
% xlabel('s_A/C_{A}');ylabel('s_R/C_{R}');
% colormap(cmap);caxis = [0 1];
% 
% subplot(2,1,2)
% h = surf(a,b,B);hold on;
% set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','linear');
% axis('square');view([0 90]);
% xlim([0 10^5]);ylim([0 10^5])
% xlabel('s_A/C_{A}');ylabel('s_R/C_{R}');
% colormap(cmap);caxis = [0 1];




%%
cmap=cbrewer('seq', 'YlOrRd', 25);
figure(2)
d = [0.1 1 100]
for i = 1:length(d)
    a = logspace(-3,3,100);
    b = logspace(-3,3,100);
    [A,B] = meshgrid(a,b);
    
    ff = (1./(1 + B./A*d(i) + B*(1+d(i))));
    
    subplot(1,3,i)
    h = contourf(a,b,ff,5);hold on;
    plot3(ones(size(a))*d(i)/(d(i)+1),b,ones(size(a))*10^5,'--k','linewidth',1);
    %         set(h,'edgecolor','none','facecolor','interp');
    set(gca,'xscale','log','yscale','log','zscale','linear');
    axis('square');view([0 90]);
    xlim([0 10^3]);ylim([0 10^3])
    %         xlabel('A/K_{A}');ylabel('R/K_{R}');
    colormap(cmap);caxis([0 1]);
    
end
Set_fig_RE(figure(2),22,22,22)

%% cooperativity

close all

a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);


aa = logspace(-2,10,10000);
bb = logspace(-2,10,10000);

na = [1 2 4];
nb = [1 2 4]
Kta = [10^10 10^2 10^2];
Ktb = [10^10 10^10 10^2];

figure(1)

for i = 1:length(na)
    for k = 1:length(nb)
        a = logspace(-2,5,100);
        b = logspace(-2,5,100);
        [A,B] = meshgrid(a,b);
        
        B = B;
        A = A;
        ff = 1./(1 + (B.^nb(k))./(A.^na(i)) +1./A.^na(i) );
        
        subplot(3,3,3*(k-1)+i)
        h = contourf(a,b,ff,5);hold on;
%         plot3(aa,1./(1-aa/d(i)),ones(size(aa))*10^5,'color',[0 0 0],'linewidth',2);
%         plot3(ones(size(a))*Kta(k),b,ones(size(a))*10^5,'--k','linewidth',1);
%         plot3(a,ones(size(a))*Ktb(k),ones(size(a))*10^5,'--k','linewidth',1);
% %         set(h,'edgecolor','none','facecolor','interp');
        set(gca,'xscale','log','yscale','log','zscale','linear');
                    set(gca,'xtick',[1  100   10000]);        set(gca,'xticklabel',[1  2   4]);
        set(gca,'ytick',[1  100   10000]);set(gca,'yticklabel',[1  2   4]);
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
%         xlabel('A/K_{A}');ylabel('R/K_{R}');
        colormap(cmap);caxis([0 1]);
    
    end
end
Set_fig_YS(figure(1),9,9,9)
%% depho - pho
% YONI: 20140804 Renan generated a color bar from here
close all
cmap=cbrewer('seq', 'YlOrRd', 25);

a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);


aa = logspace(-2,10,10000);
bb = logspace(-2,10,10000);

r = [10 1 0.1];


figure(1)

for i = 1:length(r)
        a = logspace(-2,5,100);
        b = logspace(-2,5,100);
        [A,B] = meshgrid(a,b);
   
        ff = (1./(1 + (B./A)*r(i) +  B.*(1+r(i)) ));
        
        subplot(1,3,i)
        h = contourf(a,b,ff,5);hold on;
        plot3(ones(size(a))*r(i)/(1+r(i)),b,ones(size(a))*10^5,'k','linewidth',1);
        plot3(a,ones(size(a))*1/(1 + r(i)),ones(size(a))*10^5,'k','linewidth',1);

%         set(h,'edgecolor','none','facecolor','interp');
        set(gca,'xscale','log','yscale','log','zscale','linear');
           set(gca,'xtick',[1  100   10000]);        set(gca,'xticklabel',[1  2   4]);
        set(gca,'ytick',[1  100   10000]);set(gca,'yticklabel',[1  2   4]);
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
%         xlabel('A/K_{A}');ylabel('R/K_{R}');
        colormap(cmap);
        caxis([0 1]);
        colorbar;
    
end

Set_fig_RE(figure(1),9,9,9)


end

