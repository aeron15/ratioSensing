% Gal computaintal mdoel

clear all
close all
clc


%% Parameters C = [nM] time = [min]

kf83 = 100;     x(1) = kf83;
kr83 = 1;       x(2) = kr83;
kf84 = 100;     x(3) = kf83;
kr84 = 10;     x(4) = kr84;

aG1 = 10;       x(5) = aG1;
aG3 = 10*0.8;       x(6) = aG3;
aG4 = 1;      x(7) = aG4;
a0G80 = 0.5;  x(8) = a0G80;
aG80 = 2.5 ;     x(9) = aG80;
amig = 2;       x(10) = amig;

KG1 = 10;        x(11) = KG1;
KG3 = 10;        x(12) = KG3;
KG80 = 10;       x(13) = KG80;

Kmig1 = 5;     x(14) = Kmig1;
Kmig3 = 5;      x(15) = Kmig3;
Kmig4 = 5;     x(16) = Kmig4;

n1 = 2;         x(17) = n1;
n3 = 2;         x(18) = n3;
n80 = 2;        x(19) = n80;

d = 0.004;
gG1 = d;        x(20) = gG1;
gG3 = d;        x(21) = gG3;
gG4 = d;        x(22) = gG4;
gG80 = d;       x(23) = gG80;
gG81 = d;       x(24) = gG81;
gG83 = d;       x(25) = gG83;
gG84 = d;       x(26) = gG84;
gMig = 0.004;   x(27) = gMig;

Kgal = 2^-5.5;  x(28) = Kgal;
Kglu = 2^-4;  x(29) = Kglu;

% Glu = 0
load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\mean_matrix_bfp_yfp');
agal_vec = 0;
agal_vec =  logspace(-3,0.3,50);

gal{1} = [0 2.^[-9:0.5:2]];
glu{1} = [0 2.^[-7:0.5:0]];
pcolor(log2(gal{1}),log2(glu{1}),flipud(mean_matrix));

[T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,0,0,1e7,x),[0 12*60],[0 0 0 0 0 0 0 0 0 0 0]); % Solve ODE
Y0 = Y(end,:);
    
for i = 1:length(agal_vec);
    agal = agal_vec(i);
    [T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,agal,0,1e7,x),[0 8*60],Y0); % Solve ODE
    G_vec(i,:) = Y(end,:);
end

close all

figure(1)
plot(log2(agal_vec),G_vec(:,[4 7 11 1 2 3 6]),'.-');
legend({ 'Gal4' 'Gal4-80' 'Mig1' 'Gal1' 'Gal80' 'Gal3' 'Gal3-80'})
Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')

figure(2)
plot(log2(agal_vec),mat2gray(G_vec(:,[1])),'.-');hold on
plot(log2(gal{1}),mat2gray(mean_matrix(end,:)),'r.-');

%% Double Gradiant


clear f_vec G4_vec RA_vec RE_vec Y YY

gal_res = [-9:0.5:2];
glc_res = [-7:0.5:0];

agal_vec = [0 2.^gal_res];
aglu_vec = [0 2.^glc_res];

e = 1e7
for i = 1:length(agal_vec);
    for j = 1:length(aglu_vec);
        [i,j]
        agal = agal_vec(i);
        aglu = aglu_vec(j);

         [T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,agal,aglu,e,x),[0 8*60],Y0); % Solve ODE
        
         for k = 1:11
             YY(i,j,k) = Y(end,k);
         end
 

    end
end
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
surf(log2(agal_vec),log2(aglu_vec),mat2gray(YY(:,:,1)'));hold on;view(0,90);
axis('square');
xlabel('Gal');ylabel('Glu');title('Gal 1')
colormap(cmap);
xlim([-9 1])

figure(3)
surf(log2(agal_vec),log2(aglu_vec),YY(:,:,11)');hold on;view(0,90);
axis('square');
xlabel('Gal');ylabel('Glu');title('Mig1')
colormap(cmap);

figure(4)
pcolor(log2(gal{1}),log2(glu{1}),flipud(mean_matrix));
axis('square');
xlabel('Gal');ylabel('Glu');title('Mig1')
colormap(cmap);
% Compute model slope

fit_cuttoff = [2^-1 2^-7];
mid_value = 2^-4;
cutoff=0.2;
i=1
[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(YY(:,:,4)',max(max(YY(:,:,4)')),min(min(YY(:,:,4)')),cutoff,agal_vec,aglu_vec,fit_cuttoff,mid_value);

figure(5);
plot(log2(x),log2(y),'ok','markersize',1.5);hold on;
plot(log2(x),s(log2(x)),'k');hold on;

a = 0.0076;
b = 3.48514;
c= 0;%-4579.66;
d = 0;% 7.86358E6;
plot(log2(a+aglu_vec.*b+aglu_vec.^2.*c+aglu_vec.^3.*d),log2(aglu_vec))
xlim([-9 1]);ylim([-7 0]);
