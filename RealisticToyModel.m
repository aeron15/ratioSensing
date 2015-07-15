clear all
close all
clc

M = load('mean_matrix_bfp_yfp.mat');
M = flipud(M.mean_matrix);
gal = [ 2.^[-9.5:0.5:2]];
glu = [2.^[-7.5:0.5:0]];

figure(1)
% subplot(1,2,1)

h = surf(log2(gal),log2(glu),M);
set(h,'edgecolor','k','facecolor','flat');
set(gca,'xscale','linear','yscale','linear','zscale','linear');
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 YFP');

%%
% gal = [0 logspace(-3,0.6021,24)];
% glu = [0 logspace(-2.1,0,16)];
[Gal,Glu] = meshgrid(gal,glu);
g4t = 100;
g80t = 100;
K84 = 0.25;
K83 = 0.01;
g3t0 = 7000/max(gal);
mig1t0 = 1000;
g1t=100;
K4DNA = 15;
KmigDNA = 15;
mig1 = (5 + mig1t0*Glu);

syms g80 g3t
s = solve(g80t==g80 + g3t*g80/(K83+g80) + g4t*g80/(K84+g80),g80) ;
s1 =subs(s(1),g3t,g3t0*Gal*1./(1+mig1/KmigDNA));
% s2 =subs(s(2),g3t,g3t0*Gal);
% figure(1)
% plot(gal,s1,gal,s2);set(gca,'xscale','log')
g4 = g4t*1./ (1 + real(s1)./K84 )*1./(1+mig1/(100000*KmigDNA));

%% Find g3t==g4t
g3=g3t0*Gal*1./(1+mig1/KmigDNA);
for i =1:length(glu)
    [val,ind_j(i)] = min(abs(g3(i,:)- g4t));
    ind_i(i) = i;
    
    [val,ind_j_2(i)] = min(abs(g4(i,:)-max(g4(i,:)/2)));
end
%%



e = 10^20;
g1 = 1./(1 + K4DNA./g4 + (K4DNA/KmigDNA).*(mig1./g4) + (1/e)*(mig1./KmigDNA) );

figure(2)
plot([-10 log2(gal(2:end))],g4(1,:));xlabel('Gal %');ylabel('Gal4');hold on;
plot([-10 log2(glu(2:end))],mig1(:,1),'r');
% plot([-10 log2(gal(2:end))],g1(:,1),'g')

figure(1)
subplot(1,2,2)
surf(log2(gal),log2(glu),g1);
axis('equal');view(0,90)
xlabel('Gal');ylabel('Glu');title('GAL1 YFP');

figure(3)
pcolor(log2(gal),log2(glu),g4);hold on;
plot(log2(gal(ind_j)),log2(glu),'xk');
plot(log2(gal(ind_j_2)),log2(glu),'xb');
plot(log2(g4t/g3t0*(1 + 2*(K4DNA/g4t)*(mig1t0/KmigDNA).*glu)),log2(glu),'xy')

figure(4)
plot(log2(gal),g4(1,:))

figure(5)
plot(log2(gal),real(s1));hold on;
plot(log2(gal),g80t./( 1 + g4t/K84 + (g3t0*gal)./K83 ) ,'b');
plot(log2(gal),g80t - g4t - g3t0*gal);

figure(6)
plot(log2(gal),g3t0*gal);hold on
plot(log2(gal),g4,'g')

figure(1)