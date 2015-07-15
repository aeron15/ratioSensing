%% Estimate Glc consumption



clear all
close all
clc


%%
td = 90; % doubling time
g = log(2)/(td/60); % 1/hr
u = 0.003 % %/doubling;

t = [0:0.01:30]; % hrs

od_0 = 2.3506e-05; % [od]
od_0 = 0.02;
c_0 = 2^-6; %
%%
r = 0.5;
figure(1)
subplot(2,1,1)
plot(t/(td/60),od_0*exp(g*t),'k');hold on;ylabel('OD');
title(['OD(0) = ',num2str(od_0)]);

subplot(2,1,2)
plot(t/(td/60),c_0 - u*(od_0)*(exp(g*t)-1),'k');ylabel('Glc');
Set_fig_YS(figure(1),12,12,18);
% plot(t,c_0 - u*0.05*t,'g');
% plot(t,od_0*exp(g*t)*u)
refline(0,c_0*0.9);

%%
t_inc = 8;

od_o = c_0*g*(1-r)/(u*(exp(g*t_inc)-1))

%%

prec = (c_0 - od_0*u*(exp(t(end)*g)-1))/c_0 
