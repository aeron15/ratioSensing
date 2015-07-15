%% Model the El-Sammed 2014 model


clear all
close all
clc


agal_vec = 0;
agal_vec = logspace(-3,0.3,50);

[T Y] = ode15s(@(t,y) MurrayEqsDetails(t,y,0,0,1),[0 12*60],[0 0 0 0 0 0 0 0 0 0 0]); % Solve ODE
Y0 = Y(end,:);
    
for i = 1:length(agal_vec);
    agal = agal_vec(i);
    [T Y] = ode15s(@(t,y) MurrayEqsDetails(t,y,agal,0,1),[0 12*60],Y0); % Solve ODE
    G_vec(i,:) = Y(end,:);
end


figure(1)
plot(log2(agal_vec),G_vec(:,[4 7 11 1 2 3 6]),'.-');
legend({ 'Gal4' 'Gal4-80' 'Mig1' 'Gal1' 'Gal80' 'Gal3' 'Gal3-80'})
Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')

figure(2)
plot(log2(agal_vec),G_vec(:,[4]),'.-');hold on
plot(log2(agal_vec),G_vec(:,[2])+G_vec(:,[7])+G_vec(:,[6])-0.5*(G_vec(:,[7])+G_vec(:,[4])),'k.-');
plot(log2(agal_vec),G_vec(:,[3])+G_vec(:,[6]),'r.-');

