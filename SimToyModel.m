%% Simulate the toy model and expend it

clear all
close all
clc



%%

[T Y] = ode15s(@(t,y) SimToyModelEqEff(t,y,0,0,1),[0 8*60],[1 1 1 1 1 1 1 1 1 1 1]); % Solve ODE
Y0 = Y(end,:);

gal = [ 2.^[-9.5:0.5:2]];
glu = [2.^[-7.5:0.5:0]];

agal_vec = [0 2.^[-9:0.5:2]];
aglu_vec = [0 2.^[-7:0.5:0]];
e = 1;
for i = 1:length(agal_vec);
    for j = 1:length(aglu_vec);
        [i,j]
        agal = agal_vec(i);
        aglu = aglu_vec(j);

         [T Y] = ode15s(@(t,y) SimToyModelEqEff(t,y,agal,aglu,e),[0 20*60],Y0); % Solve ODE
        
         for k = 1:11
             YY(i,j,k) = Y(end,k);
         end
 

    end
end

%% 

% y(1) = G1;
% y(2) = G80;
% y(3) = G3;
% y(4) = G4;
% y(5) = G80G1
% y(6) = G80G3
% y(7) = G80G4
% y(8) = EA (E G4)
% y(9) = RE
% y(10) = RA
% y(11) = mig

%%
close all

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
load gal_dose
load glu_dose

figure(2)
subplot(1,2,1)
h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 Log(YFP) [AU]');
colorbar
% caxis([0 1]);

figure(2)
subplot(1,2,2)
h = surf(log2(gal_dose/5.5e9),log2(glu_dose/5.5e9),Mx);hold on;
set(h,'edgecolor','k','facecolor','flat');set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view([0 90]);
xlabel('Gal');ylabel('Glu');title('GAL1 Bennet model 8 hrs [molecules/cell]');
ylim([-7 0])
colorbar
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
%%
figure(2)
plot(log2(agal_vec),YY(:,1,4)');hold on;
xlabel('Gal');ylabel('Gal4');

figure(3)
plot(log2(aglu_vec),YY(1,:,11)');hold on;
xlabel('Glu');ylabel('Mig1');

figure(4)
plot(log2(agal_vec),squeeze(YY(:,1,[2 4 7])));hold on;
legend({'Gal80','Gal4','Gal84'});
xlabel('Gal');ylabel('Gal');


figure(5)
plot(log2(agal_vec),squeeze(YY(:,1,[3 6])));hold on;
xlabel('Gal');ylabel('Gal');
legend({'Gal3','Gal83'});

figure(6)
plot(log2(agal_vec),squeeze(YY(:,1,[8 9 10])));hold on;
xlabel('Gal');ylabel('Gal');
legend({'EA','RE','RA'});
