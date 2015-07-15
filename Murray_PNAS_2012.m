%% Murray Gal Model based on Venturelli et. al PNAS 2012

clear all
close all
clc


% Effctive
close all
% for i = 1:length(agal_vec);
%         agal = agal_vec(i);
%         [T Y] = ode15s(@(t,y) MurrayEqs(t,y,agal),[0 8*60],[0 100 0 0]); % Solve ODE
%         G4_vec(i) = Y(end,1);
% end
% figure(1);hold on;
% plot(agal_vec,G4_vec,'-ok');
% Set_fig_YS(figure(1),12,12,12);

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

agal_vec = 0;
agal_vec = 0;% logspace(-3,0.3,50);

[T Y] = ode15s(@(t,y) MurrayEqsDetails(t,y,0,0,1),[0 12*60],[0 0 0 0 0 0 0 0 0 0 0]); % Solve ODE
Y0 = Y(end,:);
    
for i = 1:length(agal_vec);
    agal = agal_vec(i);
    [T Y] = ode15s(@(t,y) MurrayEqsDetails(t,y,agal,0,1),[0 12*60],Y0); % Solve ODE
    G_vec(i,:) = Y(end,:);
end

close all

figure(1)
plot(log2(agal_vec),G_vec(:,[4 7 11 1 2 3 6]),'.-');
legend({ 'Gal4' 'Gal4-80' 'Mig1' 'Gal1' 'Gal80' 'Gal3' 'Gal3-80'})
Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')

figure(2)
plot(log2(agal_vec),G_vec(:,[4]),'.-');hold on
plot(log2(agal_vec),G_vec(:,[2])+G_vec(:,[7])+G_vec(:,[6])-0.5*(G_vec(:,[7])+G_vec(:,[4])),'k.-');
plot(log2(agal_vec),G_vec(:,[3])+G_vec(:,[6]),'r.-');

% figure(1);
% plot(T,Y(:,8:10),'.-');
% legend({ 'EA' 'RE','RA'})
% Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')
% 
% 
% figure(2);
% plot(T,Y(:,[4 11]),'.-');
% legend({ 'Gal4' 'Mig1'})
% Set_fig_YS(figure(2),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')
% 
% 
% figure(3);
% plot(T,Y(:,[ 2 4 5 6 7]),'.-');
% legend({ 'Gal80' 'GAl4' 'Gal801' 'Gal803' 'Gal804'})
% Set_fig_YS(figure(3),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]');
% 
% 
% figure(4);
% plot(T,Y(:,[1 2 3 4]),'.-');
% legend({ 'Gal1' 'GAl80' 'Gal3' 'GAl4'})
% Set_fig_YS(figure(4),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')
% 
% 
% %
% legend({'1' '80' '3' '4' '81' '83' '84' 'EA' 'RE','RA' 'mig1'})
% Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')
% 
% figure(5);
% plot(agal_vec,G_vec,'.-')
% legend({'1' '80' '3' '4' '81' '83' '84' 'EA' 'RE','RA' 'mig1'})
% Set_fig_YS(figure(5),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')

%%

clear f_vec G4_vec RA_vec RE_vec Y YY
% agal_vec = [0 logspace(-3,0.3,20)];
% aglu_vec = [0 logspace(-2,0.3,20)];

gal_res = [-9:0.5:2];
glc_res = [-7:0.5:0];


agal_vec = [0 2.^gal_res];
aglu_vec = [0 2.^glc_res];

for i = 1:length(agal_vec);
    for j = 1:length(aglu_vec);
        [i,j]
        agal = agal_vec(i);
        aglu = aglu_vec(j);

         [T Y] = ode15s(@(t,y) MurrayEqsDetails(t,y,agal,aglu,1e7),[0 12*60],Y0); % Solve ODE
        
         for k = 1:11
             YY(i,j,k) = Y(end,k);
         end
 

    end
end
agal_vec(1) = 2^-10;
aglu_vec(1) = 2^-8;
    %%
close all
cmap = buildcmap('bwr');
cmap = 'jet';

figure(1)
surf(log2(agal_vec),log2(aglu_vec),YY(:,:,4)');hold on;view(0,90);
axis('square');
xlabel('Gal');ylabel('Glu');title('Gal 4');
colormap(cmap);

figure(2)
surf(log2(agal_vec),log2(aglu_vec),YY(:,:,1)');hold on;view(0,90);
axis('square');
xlabel('Gal');ylabel('Glu');title('Gal 1')
colormap(cmap);


figure(3)
surf(log2(agal_vec),log2(aglu_vec),YY(:,:,11)');hold on;view(0,90);
axis('square');
xlabel('Gal');ylabel('Glu');title('Mig1')
colormap(cmap);

%%
% 
% 
% close all
% cmap = 'jet';
% 
% figure(1)
% h = surf(log2(agal_vec),log2(aglu_vec),G1_vec);hold on;
% set(h,'edgecolor','none','facecolor','interp');
% axis('square');view([0 90]);xlim([10^-0.5 10^0.3]);ylim([10^-2 10^0.3]);
% xlabel('Gal');ylabel('Glu');
% colormap(cmap);
% 
% figure(2)
% h = surf(log2(agal_vec),log2(aglu_vec),G4_vec);hold on;
% set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
% axis('square');view([0 90]);
% xlabel('Gal');ylabel('Glu');
% colormap(cmap);
% 
% figure(3)
% h = surf(agal_vec,aglu_vec,G4_vec);hold on;
% set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
% axis('square');view([0 90]);
% xlabel('Gal');ylabel('Glu');
% colormap(cmap);
% 
% 
% load data_e_1e5
% figure(2)
% h = surf(agal_vec,aglu_vec,G1_vec);hold on;
% set(h,'edgecolor','none','facecolor','interp');set(gca,'xscale','log','yscale','log','zscale','log');
% axis('square');view([0 90]);
% xlabel('Gal');ylabel('Glu');xlim([10^-0.5 10^0.3]);ylim([10^-2 10^0.3]);
% colormap(cmap);


%%

