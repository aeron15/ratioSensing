% Fit a model to the data

close all
clear all
clc

M = load('mean_matrix_bfp_yfp.mat');
M = flipud(M.mean_matrix);
gal = [ 2.^[-9.5:0.5:2]];
glu = [2.^[-7.5:0.5:0]];
load Mx
load gal_dose
load glu_dose

d = 0.004;
amig = 1;       x0(1) = amig;
aG4 = 0.5;        x0(2) = aG4;
aG80 = 3;       x0(3) = aG80;
aG3 = 7;       x0(4) = aG3;
aG1 = 1;        x0(5) = aG1;

K4DNA = 20;     x0(6) = K4DNA;
KmigDNA = 20;   x0(7) =KmigDNA; 

kr84 = 25;       x0(8) = kr84;
kf84 = 100;     x0(9) = kf84;
kr83 = 1;       x0(10) = kr83;
kf83 = 100;     x0(11) = kf83;
kgal = -4;      x0(12) = kgal;
amig0 = 0.01;      x0(13) = amig0;

e = 1e7;

% options = optimset('TolFun',10,'TolX',1e-8);
% [F,val] = fsolve(@ModelFit,x0,options);




%

figure(100)
subplot(1,2,1)
h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 YFP');

figure(100)
subplot(1,2,2);hold on
plot(log2(gal),mat2gray(M(1,:)),'.-k');
plot(log2(gal_dose/5.5e9),mat2gray(Mx(1,:)),'.-r');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('square');
xlabel('Gal');ylabel('Gal4 YFP');

Set_fig_YS(figure(100),12,12,12)
%% Run the model

% initiail conditions
[T Y] = ode15s(@(t,y)... 
SimToyModelEqEff(t,y,0,0,e,...
                                d,x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8), x0(9), x0(10), x0(11),x0(12),x0(13)),...
                               [0 20*60],0.1*ones(1,11)); % Solve ODE
Y0 = Y(end,:);
% Glu = 0;

agal_vec = [0 2.^[-9:0.5:2]];
aglu_vec = [0 2.^[-7:0.5:0]];

for i = 1:length(agal_vec);
    for j = 1:length(aglu_vec);
        [i,j]
        agal = agal_vec(i);
        aglu = aglu_vec(j);

         [T Y] = ode15s(@(t,y)... 
                 SimToyModelEqEff(t,y,agal,aglu,e,...
                                d,x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8), x0(9), x0(10), x0(11),x0(12),x0(13)),...
                               [0 10*60],Y0); % Solve ODE;
        
         for k = 1:11
             YY(i,j,k) = Y(end,k);
         end
 

    end
end

plot(log2(agal_vec),mat2gray(squeeze(YY(:,1,1))),'.-g');
figure
plot(log2(agal_vec),(squeeze(YY(:,1,4))),'.-g');

%%



M = load('mean_matrix_bfp_yfp.mat');
M = flipud(M.mean_matrix);
gal = [ 2.^[-9.5:0.5:2]];
glu = [2.^[-7.5:0.5:0]];

figure(1)
subplot(1,2,1)
h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 YFP');
% caxis([0 1]);

figure(1)
subplot(1,2,2)
h = surf(log2(gal),log2(glu),(YY(:,:,1)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title(['GAL1 model \epsilon = ',num2str(e)]);
% caxis([0 1]);
Set_fig_YS(figure(1),12,12,12)

load Mx
Mx = Mx;
load Mx_II
Mx_II = Mx;
load gal_dose
load glu_dose

figure(2)
subplot(1,3,1)
h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 Log(YFP) [AU]');
% colorbar
% caxis([0 1]);

figure(2)
subplot(1,3,2)
h = surf(log2(gal_dose/5.5e9),log2(glu_dose/5.5e9),Mx);hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('GAL1 Bennet model 8 hrs [molecules/cell]');
% colorbar
% caxis([0 1]);

figure(2)
subplot(1,3,3)
h = surf(log2(gal),log2(glu),(YY(:,:,1)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('GAL1 Model YS [molecules/cell]');
% colorbar
% caxis([0 1]);

Set_fig_YS(figure(2),12,12,12)


figure(3)
subplot(1,2,1)
h = surf(log2(gal),log2(glu),(YY(:,:,11)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Mig1 [molecules/cell]');
colorbar

subplot(1,2,2)
h = surf(log2(gal),log2(glu),(YY(:,:,4)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Gal4 [molecules/cell]');
colorbar

Set_fig_YS(figure(3),12,12,12)

figure(4)
subplot(2,2,1)
h = surf(log2(gal),log2(glu),(YY(:,:,4)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Gal4 [molecules/cell]');
colorbar
subplot(2,2,2)
h = surf(log2(gal),log2(glu),(YY(:,:,2)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Gal80 [molecules/cell]');
colorbar
subplot(2,2,3)
h = surf(log2(gal),log2(glu),(YY(:,:,3)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Gal3 [molecules/cell]');
colorbar
subplot(2,2,4)
h = surf(log2(gal),log2(glu),(YY(:,:,1)'));hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('Gal1 [molecules/cell]');
colorbar
%%
figure(5)
plot(log2(agal_vec),YY(:,1,4)');hold on;
xlabel('Gal');ylabel('Gal4');

figure(6)
plot(log2(aglu_vec),YY(1,:,11)');hold on;
xlabel('Glu');ylabel('Mig1');

figure(7)
plot(log2(agal_vec),squeeze(YY(:,1,[2 4 7])));hold on;
legend({'Gal80','Gal4','Gal84'});
xlabel('Gal');ylabel('Gal');


figure(8)
plot(log2(agal_vec),squeeze(YY(:,1,[3 6])));hold on;
xlabel('Gal');ylabel('Gal');
legend({'Gal3','Gal83'});