%% Parse data

clear all
close all
clc


load('C:\Users\ys151\Desktop\WorkingZone\Gal\LowODMeta\gal8hrs');

figure(1)
eu = squeeze(C(2,1,:))';
ed = squeeze(C(1,1,:))';
c = [log2(2.^[-10:1])];
errorbar(c,m,m-ed,eu-m,'ro');
gal_g = m;

load('C:\Users\ys151\Desktop\WorkingZone\Gal\LowODMeta\glc8hrs');

figure(1);hold on
eu = squeeze(C(2,1,:))';
ed = squeeze(C(1,1,:))';
c = [log2(2.^[-10:1])];
errorbar(c,m,m-ed,eu-m,'ko');
glc_g = m;



%%
xx = logspace(-3,0.3,100);
yy_gal = interp1(2.^c,gal_g,xx);
yy_glc = interp1(2.^c,glc_g,xx);

figure(2)
plot(log2(xx),yy_gal,'r');hold on;
plot(log2(xx),yy_glc,'k')


%%
cmap =  buildcmap_YS([0 0 0.5;1 1 1;0.5 0 0]);
colormap(cmap);
[s1,s2] = meshgrid(xx,xx);

for i = 1:length(xx)
    for j = 1:length(xx)
        
        if yy_glc(i)>=yy_gal(j)
            F(i,j) = yy_glc(i);
            ind(i,j) = 0;
        else
            F(i,j) = yy_gal(j);
            ind(i,j) = 1;
        end
    end
end


figure(3)
subplot(1,2,2)
h = surf(log2(xx),log2(xx),F);hold on;
set(h,'edgecolor','none','facecolor','interp');
view(0,90);axis square
xlim([-8 1]);ylim([-6 0]);
cmap =  buildcmap_YS([0 0 0.5;1 1 1;0.5 0 0]);

subplot(1,2,1)
h = surf(log2(xx),log2(xx),ind);hold on;
set(h,'edgecolor','none','facecolor','interp');
view(0,90);axis square
xlim([-10 1]);ylim([-10 1]);
colormap(cmap)
