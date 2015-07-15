% Experimental fitness experiment
clear all
close all
clc

color = [0 0 0;...
    0.85         0.16    0;...
    0 0.5 0;...
    0    153/255    1];

M = load('mean_matrix_bfp_yfp.mat');
M = flipud(M.mean_matrix);
gal = [ 2.^[-9.5:0.5:2]];
glu = [2.^[-7.5:0.5:0]];

figure(1)
h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 YFP');

%%
load ('C:\Users\ys151\Dropbox\20131120\output\plates_bfp.mat')
load ('C:\Users\ys151\Dropbox\20131120\output\plates_mCh.mat')
close all

media = {'2^{-4.5}%Glc 0% Gal' '2^{-6}%Glc 0% Gal' '2^{-4.5}%Glc 2% Gal' '2^{-6}%Glc 2% Gal'};
t = [15 17 19 21];
% No gal

x1 =   struct2cell(plates_mCh.Flask_1);
x2 =   struct2cell(plates_mCh.Flask_2);
x3 =   struct2cell(plates_mCh.Flask_3);
x4 =   struct2cell(plates_mCh.Flask_4);
x5 =   struct2cell(plates_mCh.Flask_5);
x6 =   struct2cell(plates_mCh.Flask_6);
x7 =   struct2cell(plates_mCh.Flask_7);
x8 =   struct2cell(plates_mCh.Flask_8);

for i =1:4
    n_wt_h(1,i) = x1{i}.Fraction/x1{1}.Fraction;
    n_wt_h(2,i) = x2{i}.Fraction/x2{1}.Fraction;
    n_wt_l(1,i) = x5{i}.Fraction/x5{1}.Fraction;
    n_wt_l(2,i) = x6{i}.Fraction/x6{1}.Fraction;
    n_wt_h_gal(1,i) = x3{i}.Fraction/x3{1}.Fraction;
    n_wt_h_gal(2,i) = x4{i}.Fraction/x4{1}.Fraction;
    n_wt_l_gal(1,i) = x7{i}.Fraction/x7{1}.Fraction;
    n_wt_l_gal(2,i) = x8{i}.Fraction/x8{1}.Fraction;
    
    wt_h(1,i) = x1{i}.Fraction;
    wt_h(2,i) = x2{i}.Fraction;
    wt_l(1,i) = x5{i}.Fraction;
    wt_l(2,i) = x6{i}.Fraction;
    wt_h_gal(1,i) = x3{i}.Fraction;
    wt_h_gal(2,i) = x4{i}.Fraction;
    wt_l_gal(1,i) = x7{i}.Fraction;
    wt_l_gal(2,i) = x8{i}.Fraction;
end

%%
figure(2);
subplot(1,2,1);hold on;title('Fraction of WT');

for i =1:2
plot(t,wt_h(i,:),'.-','color',color(1,:))
plot(t,wt_l(i,:),'.-','color',color(2,:))
plot(t,wt_h_gal(i,:),'.-','color',color(3,:))
plot(t,wt_l_gal(i,:),'.-','color',color(4,:))
xlabel('time [hr]');
if i==1
legend(media)
end
end
set(gca,'YGrid','on');

subplot(1,2,2);hold on
title('Normalized fraction of WT');
plot(t,log(n_wt_h),'.-','color',color(1,:))
plot(t,log(n_wt_l),'.-','color',color(2,:))
plot(t,log(n_wt_h_gal),'.-','color',color(3,:))
plot(t,log(n_wt_l_gal),'.-','color',color(4,:))
xlabel('time [hr]');ylabel('log(r(t)/r(0))')
set(gca,'YGrid','on');

Set_fig_YS(figure(2),12,12,12);

figure(3);hold on;
plot(t,log(n_wt_l_gal),'.-')
s1 = fit(t',log(n_wt_l_gal(1,:))','a*(x-15)');
s2 = fit(t',log(n_wt_l_gal(2,:))','a*(x-15)');
plot(t,s1(t),t,s2(t));
xlabel('time [hr]');ylabel('log(r(t)/r(0))')
set(gca,'YGrid','on');
title('\Deltag = 0.03-0.04 1/hr')
Set_fig_YS(figure(3),12,12,12);

