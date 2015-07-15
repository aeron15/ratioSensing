close all
clear all
clc

dir_name = 'C:\Users\ys151\Dropbox\20130628\data\'
load([dir_name,'Mx_A10_130628']);
A10 = Mx;
load([dir_name,'Mx_H3_130628']);
H3 = Mx;
H3(H3==0)=2;
%% Plot the data
close all

c_gal = [0 2.^(-9:0.5:2)]
c_glu =[0 2.^(-10:1)];
L = max( [max(A10(:)),max(H3(:))]);
figure(1)
subplot(3,2,1);
h = pcolor(c_gal,c_glu,A10);
set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
axis('square');
xlabel('Gal %');ylabel('Glu %');
caxis([2 L]);colorbar;title('A10')
Set_fig_YS(figure(1),12,12,12)

subplot(3,2,2);
h = pcolor(c_gal,c_glu,H3);
set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
axis('square');
xlabel('Gal %');ylabel('Glu %');
caxis([2 L]);colorbar;title('H3')
Set_fig_YS(figure(1),12,12,12);

[Gal,Glu] = meshgrid(c_gal,c_glu);
r = Gal./Glu;

subplot(3,2,3);
h = pcolor(r,Gal,A10);
set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
axis('square');
xlabel('Gal/Glu');ylabel('Glu %');
caxis([2 L]);colorbar
Set_fig_YS(figure(1),12,12,12)

subplot(3,2,4);
h = pcolor(r,Gal,H3);
set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
axis('square');
xlabel('Gal/Glu');ylabel('Glu %');
caxis([2 L]);colorbar;
Set_fig_YS(figure(1),12,12,12);

figure(2)
subplot(1,2,1);
h = plot(r,A10,'.-');
set(gca,'xscale','log','yscale','linear');
axis('square');
xlabel('Gal/Glu');ylabel('Glu %');
caxis([2 L]);
Set_fig_YS(figure(2),12,12,12)

subplot(1,2,2);
h = plot(r,H3,'.-');
set(gca,'xscale','log','yscale','linear');
axis('square');
xlabel('Gal/Glu');ylabel('Glu %');
caxis([2 L]);
Set_fig_YS(figure(2),12,12,12);

% %% Fit the data to the ratio
% for i = 1:length(c_gal);
%     s = fit(r(2:end,i),H3(2:end,i),'a/(1/x+b)');
%     t_H3(i) = s.b;
%     
%     s = fit(r(2:end,i),A10(2:end,i),'a/(1/x+b)');
%     t_A10(i) = s.b;
%    
%    
% end
%%


%%





%%
close all
cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
r = logspace(log(min(a)/max(b)),log(max(a)/min(b)),100);

[A,B] = meshgrid(a,a);
[A,R] = meshgrid(a,r);
% One binding site
f = log(1./(1 + B./A +1./A ));
f1 = log( 1./(1 + 1./R + 1./A ) );


c_gal = [0 2.^(-9:0.5:2)]
c_glu =[0 2.^(-10:1)];



figure(1)
h = surf(a,1./b,f);hold on;
set(h,'edgecolor','none','facecolor','interp');
plot3(a,1-a,ones(size(a))*10^5,'--k');
set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}');ylabel('B/K_{B}');
colormap(cmap);

figure(2)
h = surf(r,a,f1');
set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 (max(r))]);ylim([0 (max(a))])
xlabel('A/K_{A}/B/K_{B}');ylabel('A/K_{A}');

Set_fig_YS(figure(1),12,18,18);
Set_fig_YS(figure(2),12,18,18);

colormap(cmap);

a = logspace(-2,5,100);
b = logspace(-2,5,100);
r = logspace(log(min(a)/max(b)),log(max(a)/min(b)),100);

[A,B] = meshgrid(a,b);
[A,R] = meshgrid(a,r);
d=1;
ff = log(1./(1 + B./A +1./A + B/d));
ff1 = log( 1./(1 + 1./R + 1./A + A./R/d ) );


figure(3)
h = surf(a,b,ff);hold on
set(h,'edgecolor','none','facecolor','interp');
plot3(a,ones(size(a)),ones(size(a))*10^5,'--k');
plot3(ones(size(a)),b,ones(size(a))*10^5,'--k');

set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}');ylabel('B/K_{B}');
colormap(cmap);


figure(4)
h = surf(r,a,ff1');
set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}/B/K_{B}');ylabel('A/K_{A}');
colormap(cmap);

Set_fig_YS(figure(3),12,18,18);
Set_fig_YS(figure(4),12,18,18);

d=10;
ff = log(1./(1 + B./A +1./A + B/d));
ff1 = log( 1./(1 + 1./R + 1./A + A./R/d ) );


figure(5)
h = surf(a,b,ff);hold on;
plot3(a,(a-1)*d./(a-d),ones(size(a))*10^5,'--k');
plot3(d*(1-b)./(d-b),b,ones(size(a))*10^5,'--k');
set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}');ylabel('B/K_{B}');
colormap(cmap);


figure(6)
h = surf(r,a,ff1');
set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}/B/K_{B}');ylabel('A/K_{A}');
colormap(cmap);

Set_fig_YS(figure(5),12,18,18);
Set_fig_YS(figure(6),12,18,18);

%%

fff = log(-((A + B + A .*B - sqrt(4 *A.*(1 + B).^2 +(A - B + A.* B).^2))./(2 *(1 + B))));



figure(5)
h = surf(a,b,fff);hold on
set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
axis('square');view([0 90]);
xlim([0 10^5]);ylim([0 10^5])
xlabel('A/K_{A}');ylabel('B/K_{B}');



%%
close all

a = logspace(-2,5,100);
b = logspace(-2,5,100);
r = logspace(log(min(a)/max(b)),log(max(a)/min(b)),100);

[A,B] = meshgrid(a,b);
[A,R] = meshgrid(a,r);

d = [1 1.5 10 10^10];
for i = 1:length(d)
    
    f = (1./(1 + B./A +1./A + B/d(i)));

    figure(1);subplot(1,4,i)
    h = surf(a,1./b,f);hold on
    set(h,'edgecolor','none');set(gca,'xscale','log','yscale','log');
    axis('square');view([0 90]);
    xlim([0 10^5]);ylim([0 10^5])
    xlabel('A/K_{A}');ylabel('B/K_{B}');

end
