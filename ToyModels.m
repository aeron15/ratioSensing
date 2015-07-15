%% Gal Toy models

clear all
close all
clc


%% Pure ratio
close all
cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);

    ff = (1./(1 + (B.^1./A)^2));
        h = surf(a,b,ff);hold on;
    set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');

  axis('square');view([0 90]);
    xlim([0 10^5]);ylim([0 10^5])
    xlabel('A');ylabel('R');
    colormap(cmap);caxis = [0 1];
    colorbar
    Set_fig_YS(figure(1),12,18,18)

%% Exclusion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clcclc%%

cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);

d = [1 10^3 10^7];
figure(1)
for i = 1:length(d)
    
    ff = (1./(1 + B./A +1./A + B/d(i)));
    subplot(1,length(d),i)
    h = surf(a,b,ff);hold on;
     plot3(a,(a-1)*d(i)./(a-d(i)),ones(size(a))*10^5,'--y','linewidth',2);
    plot3(d(i)*(1-b)./(d(i)-b),b,ones(size(a))*10^5,'--y','linewidth',2);
     plot3(a,1./(1-a/d(i)),ones(size(a))*10^5,'--y','linewidth',2);

    set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','linear');
    axis('square');view([0 90]);
    xlim([0 10^5]);ylim([0 10^5])
    xlabel('A/K_{A}');ylabel('R/K_{R}');
    colormap(cmap);caxis = [0 1];
end
Set_fig_YS(figure(1),12,18,18)
%% Transport and deplition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure(2)

cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);

[A,B] = meshgrid(a,b);

d = [10^10];
Kta = [10^2 10^2 10^4];
Ktb = [10^2 10^4 10^10];

for k= 1:length(Kta)
    for i = 1:length(d)
        [A,B] = meshgrid(a,b);
        
        B = B*Ktb(k)./(Ktb(k) + B);
        A = A*Kta(k)./(Kta(k) + A);
        ff = (1./(1 + B./A +1./A + B/d(i)));
        subplot(length(d),length(Kta),(k-1)*length(d)+i)
        h = surf(a,b,ff);hold on;
        plot3(a,(a-1)*d(i)./(a-d(i)),ones(size(a))*10^5,'--y','linewidth',2);
        plot3(d(i)*(1-b)./(d(i)-b),b,ones(size(a))*10^5,'--y','linewidth',2);
        
        plot3(ones(size(a))*Kta(k),b,ones(size(a))*10^5,'--k');
        plot3(a,ones(size(a))*Ktb(k),ones(size(a))*10^5,'--k');
        
        set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
        xlabel('A/K^{*}_{A}');ylabel('R/K^{*}_{R}');
        colormap(cmap);caxis = [0 1];
        
    end
end

Set_fig_YS(figure(2),12,18,18);

%% feedback
close all
cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,3,100);
b = logspace(-2,3,100);
[A,B] = meshgrid(a,b);
syms AE alpha beta x d gf

%
% y1 = ( subs(s(2),{'alpha','beta','d'},{a,0,1}) );
% y2 = ( subs(s(2),{'alpha','beta','d'},{a,0,d(i)}) );
% figure(1);hold all;plot(a,y1);set(gca,'xscale','log','yscale','log');

gf = [0 100];
gf_m = max(gf);

for k = 1:length(gf)
    

    
    
    x = alpha * ( (1 + gf(k)*AE)/(1 + gf(k)) ) ;
    
    s = solve(AE == 1/( 1 + 1/x + beta/d + beta/x),AE);
    
    d_vec = [1 10 10^7];
    
    for i = 1:length(d_vec)
        
        if k==1
        ff1 =  ( subs(s(1),{'alpha','beta','d'},{A,B,d_vec(i)}) );
        else
            
        ff1 =  ( subs(s(2),{'alpha','beta','d'},{A,B,d_vec(i)}) );
        end
        ff1(find(imag(ff1)))=0;
        figure(1)
        subplot(length(gf),length(d_vec),(k-1)*length(d_vec)+i)
        h = surf(a,b,ff1);hold on;
        plot3(a,(a-1)*d_vec(i)./(a-d_vec(i)),ones(size(a))*10^5,'--y','linewidth',2);
        plot3(d_vec(i)*(1-b)./(d_vec(i)-b),b,ones(size(b))*10^5,'--y','linewidth',2);
        
        set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
        axis('square');view([0 90]);
        xlim([0 10^3]);ylim([0 10^3])
        xlabel('A/K_{A}');ylabel('R/K_{R}');
        colormap(cmap);caxis = [0 1];
        figure(2);hold all;
        plot(a,(ff1(1,:)),'.-');set(gca,'xscale','log','yscale','linear');
    end
    
end

figure(2);set(gca,'xscale','log');xlabel('A/K_{A}');ylabel('Response');
title('R=0')
Set_fig_YS(figure(2),12,18,18);

%% All together
close all
%cmap = redbluecmap;%
buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);
syms AE alpha beta x d gf

d_vec = [10^10]; % high exculsion
i = 1;



Kta = [10^10]; 
Ktb = [10^10];
j = 1;



gf = [10 10]
gf_m = 10;
k = 2;

m=1;
h=1;
vb = Kta/h;
va = Kta/m;
z=1;

    x = alpha * ( (1 + gf(k)*AE)/(1 + gf_m) ) ;
    B = B*vb./(Ktb(j) + B);
    A = A*va./(Kta(j) + A);

    s = solve(AE == 1/( 1 + 1/x + beta/d + beta/x),AE);
    
    
    for i = 1:length(d_vec)
        
        if k==1
        ff1 =  ( subs(s(1),{'alpha','beta','d'},{A,B,d_vec(i)}) );
        else
            
        ff1 =  ( subs(s(2),{'alpha','beta','d'},{A,B,d_vec(i)}) );
        end
        ff1(find(imag(ff1)))=0;
        m = zeros(size(ff1));m1 = m;m2 = m;
        ind_1 = find(ff1>0.1);ind_2 = find(ff1<0.9);
        m1(ind_1) = 1;m2(ind_2) = 1;m = ff1;
        
        figure(z)
        h = surf(a,b,ff1);hold on;
%         plot3(a,(a-1)*d_vec(i)./(a-d_vec(i)),ones(size(a))*10^5,'--y','linewidth',2);
%         plot3(d_vec(i)*(1-b)./(d_vec(i)-b),b,ones(size(b))*10^5,'--y','linewidth',2);
%          plot3(ones(size(a))*Kta(j),b,ones(size(a))*10^5,'--k');
%         plot3(a,ones(size(a))*Ktb(j),ones(size(a))*10^5,'--k');
        set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
        axis('square');view([0 90]);
        xlim([0 10^5]);ylim([0 10^5])
        xlabel('A');ylabel('R');
        colormap(cmap);caxis = [0 1];grid off;

%         figure(2);hold all;
%         plot(a,(ff1(1,:)),'.-');
    end
                  

Set_fig_YS(figure(z),12,18,18);



%% Fitness model ti = 0
clear all
close all
g = 3;
dg = 0.1;
dd = 1;
d = 0;

 cmap = buildcmap('bwr');
%cmap = redbluecmap;
n = 200;
s1 = logspace(-1,1.5,n);
s2 = logspace(-1,1.5,n);
u = 2;
t = 1;

[S1,S2] = meshgrid(s1,s2);


figure(10);hold all;
plot(s2,s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)) );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
plot(s2,u*t*s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'r' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');

%
% set(gca,'xscale','log','yscale','log');
Set_fig_YS(figure(10),12,12,12);
%


fs1 = (g+dg)./(1+1./S2) - d;
fs1(find(fs1<0))=0;

fs2 = g./(1+1./S1) - d - dd;
fs2(find(fs2<0))=0;

figure(11);hold on;
plot(s1,fs1(:,1),'k');
plot(s2,fs2(1,:),'r');

%%

% i = s1, j = s2,

for i = 1:n
    for j = 1:n
        
        [f(i,j),ind(i,j)] = max([fs1(i,j),fs2(i,j)]);

        if s1(i) >= u*t * s2(j).*(g-dd*(1+1./s2(j)))./(g+(dg+dd)*(1+s2(j)))
            f2(i,j) = fs1(i,j);
            ind2(i,j) = 1;
        else
             f2(i,j) = fs2(i,j);
            ind2(i,j) = 2;
        end
        
    end
end
%
close all
figure(1);
subplot(1,2,1)
% h=pcolor(s2,s1,f);colormap(cmap);
h = contourf(s2,s1,f,20,'edgecolor','none');colormap(cmap);;title('Fitness');

% set(h,'edgecolor','none');
axis square;hold on;


plot(s2,s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');

figure(1);
subplot(1,2,2);
h=pcolor(s2,s1,ind);colormap(cmap);title('Induction');
set(h,'edgecolor','none');axis square;hold on;
plot(s2,s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');

Set_fig_YS(figure(1),12,18,12);

% t = 1;
figure(2);
subplot(1,2,1)
% h=pcolor(s2,s1,f);colormap(cmap);
h = contourf(s2,s1,f2,20,'edgecolor','none');colormap(cmap);;title('Fitness');

% set(h,'edgecolor','none');
axis square;hold on;

plot(s2,+s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');
plot(s2,u*t*s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');

figure(2);
subplot(1,2,2);
h=pcolor(s2,s1,ind2);colormap(cmap);title('Induction');
set(h,'edgecolor','none');axis square;hold on;
plot(s2,+s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');
plot(s2,u*t*s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2)),'k' );xlabel('\sigma_{2}');ylabel('\sigma_{1}');
set(gca,'xscale','log','yscale','log');

Set_fig_YS(figure(2),12,18,12);

%% Analyze growth
close all
load Growth_Rate % A10 A10 H3 H3
load final_yield_plate_1_A10

glc = 2.^[-5:1];
glc = [0 glc]
gal = 2.^[-8:1];
gal = [0 gal]

figure(1)
subplot(2,1,1);hold on;
plot(glc,fs1_e(1:2,:),'.-k');title('A10')
plot(gal,fs2_e(1:2,:),'.-r');
set(gca,'xscale','log')
subplot(2,1,2);hold on;
plot(glc,fs1_e(3:4,:),'.-k');title('H3')
plot(gal,fs2_e(3:4,:),'.-r');
set(gca,'xscale','log')

%% Based on A10 Yeild

figure(1);hold on;

plot(glc,yield(:,1),'.k'); 
plot(gal,yield(1,1:11),'.r')
figure(2)
pcolor(yield)
%%


for i =1:4
    fs1_e(i,:) = growth_rate(1:8,1,i);
    fs2_e(i,:) = growth_rate(1,1:11,i);

end

figure(1)
subplot(1,2,1)
plot(glc,fs1_e(3:4,:),'.-');
subplot(1,2,2)
plot(gal,fs2_e(3:4,:),'.-');

%H3

figure(2)
subplot(2,1,1)
plot(glc,fs1_e(3:4,:),'.-k');hold on;

subplot(2,1,2)
plot(gal,fs2_e(3:4,:),'.-k');hold on;

%%
close all

%H3
cutoff = 4;
c20 = gal(4);
for i = 4
    
    s_glc = fit(glc',fs1_e(i,1:end)','g*x/(x+K1)');
    
    figure(i);
    subplot(1,2,1)
    hold on;
    plot(glc,fs1_e(i,:),'.-k');plot(glc,s_glc(glc),'r');
% set(gca,'xscale','log','yscale','log');
y = fs2_e(i,4:end)';x = gal(4:end)';
y(1) = 0;
s0 = x(1);
s_gal = fit(x,y,'g*(x-0.0156)/((x-0.0156)+K1)');
    
    subplot(1,2,2)
    hold on;
    plot(x,y,'.-k');plot(x,s_gal(x),'r');
%     set(gca,'xscale','log','yscale','log');

end
%%
close all
s1 = [0 logspace(-11,1,500)]
s2 = [0 logspace(-11,1,500)]

[S1,S2] = meshgrid(s1,s2);

for i = 1:length(s1)
    for j = 1:length(s2)
        
        fs1(i,j) = 0.5947*s1(i)/(s1(i) + 0.0087);
        if s2(j) <= 0.0156
            fs2(i,j) = 0;
        else
            fs2(i,j) = 0.35*(s2(j)-0.0156)/(s2(j)-0.0156 + 0.0041);
        end
        [f_e(i,j),ind_e(i,j)] = max([fs1(i,j),fs2(i,j)]);
    end
end
[u,v]=find(diff(ind_e,1,2));
v = v+1;
figure;
h=pcolor(log2(s2),log2(s1),ind_e);title('Induction');
set(h,'edgecolor','none');axis square;hold on;
plot(s2(v),s1(u),'k')
% set(gca,'xscale','log','yscale','log');

figure;
h=pcolor(s2,s1,fs1);title('Induction');
set(h,'edgecolor','none');axis square;hold on;
set(gca,'xscale','log','yscale','log');

figure;
h=pcolor(s2,s1,fs2);title('Induction');
set(h,'edgecolor','none');axis square;hold on;
set(gca,'xscale','log','yscale','log');

figure(2)
subplot(2,1,1)
plot(log2(glc),fs1_e(3:4,:),'.-k');hold on;
plot(log2(s1),fs1(1:end,1))
subplot(2,1,2)
plot(log2(gal),fs2_e(3:4,:),'.-k');hold on;
plot(log2(s2),fs2(1,1:end))
figure
hold on 
plot(log2(s1),fs1(1:end,1))
plot(log2(s2),fs2(1,1:end))


%%    

g = 0.3;
dg = 0.3;
dd = 0.00;
d = 0;
K1 = 0.04;
K2 = 0.056;

% s1 = linspace(0.1,0.5,n);
s2 = logspace(-15,1,50)


figure(10);hold all;
plot(log2(s2),log2(s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2))) ,'.-k');xlabel('\sigma_{2}');ylabel('\sigma_{1}');
% xlim([-8 2]);ylim([-8 2]);
% set(gca,'xscale','log','yscale','log');

%% Plot the real induction  data
dir_name = 'C:\Users\ys151\Dropbox\20130628\data\'
load([dir_name,'Mx_A10_130628']);
A10 = Mx;
load([dir_name,'Mx_H3_130628']);
H3 = Mx;
H3(H3==0)=2;

close all

c_gal = [-9.5 (-9:0.5:2)];
c_glu =[-11 (-10:1)];

L = max( [max(A10(:)),max(H3(:))]);
% A10


figure(1)
h = pcolor(c_gal,c_glu,A10);
set(h,'edgecolor','none');
axis('square');
xlabel('Gal %');ylabel('Glu %');
caxis([2 L]);colorbar;title('A10');hold on;
% plot(log2(s2),log2(s2.*(g-dd*(1+1./s2))./(g+(dg+dd)*(1+s2))) ,'.-k');
% set(gca,'xtick',[1:24],'xticklabel',c_gal,'ytick',[1:13],'yticklabel',c_glu);

Set_fig_YS(figure(1),12,12,12)

% H3

figure(2)
h = pcolor(c_gal,c_glu,H3);
set(h,'edgecolor','none');
% set(gca,'xtick',[1:24],'xticklabel',c_gal,'ytick',[1:13],'yticklabel',c_glu);

% set(gca,'xscale','log','yscale','log');
axis('square');
xlabel('Gal %');ylabel('Glu %');
caxis([2 L]);colorbar;title('H3');hold on;

plot(log2(s2(v)),log2(s1(u)),'.-k');

Set_fig_YS(figure(2),12,12,12);
%% Direct reprssion of A

cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[A,B] = meshgrid(a,b);

d = [1 10 10^7];
figure(1)
for i = 1:length(d)
    A = A./(1+B);
    ff = (1./(1 + B./A +1./A + B/d(i)));
    subplot(1,length(d),i)
    h = surf(a,b,ff);hold on;
     plot3(a,(a-1)*d(i)./(a-d(i)),ones(size(a))*10^5,'--y','linewidth',2);
    plot3(d(i)*(1-b)./(d(i)-b),b,ones(size(a))*10^5,'--y','linewidth',2);
    set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','linear');
    axis('square');view([0 90]);
    xlim([0 10^5]);ylim([0 10^5])
    xlabel('A/K_{A}');ylabel('R/K_{R}');
    colormap(cmap);caxis = [0 1];
end
Set_fig_YS(figure(1),12,18,18)
%% competition on the trasporter level


cmap = redbluecmap;%   buildcmap('bwr');
a = logspace(-2,5,100);
b = logspace(-2,5,100);
[AA,BB] = meshgrid(a,b);

d = [1 10 10^7];
figure(1)

for i = 1:length(d)
    A = 1./(1 + 1./AA +BB./AA) ;
    B = 1./(1 + 1./BB +AA./BB) ;

    
    ff = (1./(1 + B./A +1./A + B/d(i)));
    subplot(1,length(d),i)
    h = surf(a,b,ff);hold on;
     plot3(a,-a*(1-1/(d(i))),ones(size(a))*10^5,'--y','linewidth',2);
    set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','linear');
    axis('square');view([0 90]);
    xlim([0 10^5]);ylim([0 10^5])
    xlabel('A/K_{A}');ylabel('R/K_{R}');
    colormap(cmap);caxis = [0 1];
end
Set_fig_YS(figure(1),12,18,18)
