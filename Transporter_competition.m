function Transporter_competition

%% Transporter uptakes.

doublings=2;

Khxt_glc = 1
Khxt_gal = 200;
vmax_hxt = 1;


Kgal2_glc = 2;
Kgal2_gal = 1;
vmax_gal2= 0;

c_max = 1;
gal_res = [-9:0.1:1];
glc_res = [-7:0.1:0];

agal_vec = [ 2.^gal_res]*55;
aglu_vec = [ 2.^glc_res]*55;

[glc_out,gal_out] = meshgrid(aglu_vec,agal_vec)

glc_in_hxt = vmax_hxt./(1 + Khxt_glc./glc_out + Khxt_glc/Khxt_gal*gal_out./glc_out);
gal_in_hxt = vmax_hxt./(1 + Khxt_gal./gal_out + Khxt_gal/Khxt_glc*glc_out./gal_out);

glc_in_gal2 = vmax_gal2./(1 + Kgal2_glc./glc_out + Kgal2_glc/Kgal2_gal*gal_out./glc_out);
gal_in_gal2 = vmax_gal2./(1 + Kgal2_gal./gal_out + Kgal2_gal/Kgal2_glc*glc_out./gal_out);

figure(1)
subplot(3,2,1)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),glc_in_hxt');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('HXT glc');caxis([0 c_max]);

subplot(3,2,2)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),gal_in_hxt');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('HXT gal');caxis([0 c_max]);

subplot(3,2,3)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),glc_in_gal2');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('GAL2 glc');caxis([0 c_max]);

subplot(3,2,4)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),gal_in_gal2');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('GAL2 gal');caxis([0 c_max]);

subplot(3,2,5)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),glc_in_gal2'+glc_in_hxt');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('Total glc');caxis([0 c_max]);

subplot(3,2,6)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),gal_in_gal2'+gal_in_hxt');set(h,'edgecolor','none');
xlabel('Gal');ylabel('Glc');title('Total gal');caxis([0 c_max]);

figure(2)
subplot(3,1,1)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),gal_in_hxt'./glc_in_hxt');set(h,'edgecolor','none');
view(0,90);xlabel('Gal');ylabel('Glc');title('HXT gal/glc');caxis([0 c_max]);

subplot(3,1,2)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),gal_in_gal2'./glc_in_gal2');set(h,'edgecolor','none');
view(0,90);xlabel('Gal');ylabel('Glc');title('Gal2 gal/glc');caxis([0 c_max]);

subplot(3,1,3)
h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),(gal_in_hxt'+gal_in_gal2')./(glc_in_hxt'+glc_in_gal2'));set(h,'edgecolor','none');
view(0,90);xlabel('Gal');ylabel('Glc');title('GAL2 glc');caxis([0 c_max]);


%% Draw the c12/c13 uptake experiment
figure(3)

R = (glc_in_hxt'+glc_in_gal2')./((gal_in_hxt'+gal_in_gal2')+(glc_in_hxt'+glc_in_gal2')); % ratio of c13/c12
% R = 1/2^doublings+(2^doublings-1)*R/2^doublings; % one doubling correct for a diffrent growth.
R = R/2;
f = 1./(1 + R);

%h = pcolor(log2(agal_vec/55),log2(aglu_vec/55),R);set(h,'edgecolor','none');
h = surf(log2(agal_vec/55),log2(aglu_vec/55),R);set(h,'edgecolor','none');

view(0,90);xlabel('Gal');ylabel('Glc');title('GAL2 glc');caxis([0 c_max]);

% S glucose -6.5 and galactose -5 ( 3 tubes)
% S glucose -6.5 and galactose -1 ( 3 tubes)    
% S glucose -6.5 and galactose 1 ( 3 tubes)
  
% S glucose -4 and galactose -5 ( 3 tubes)
% S glucose -4 and galactose -1 ( 3 tubes)
% S glucose -4 and galactose 1 ( 3 tubes)
glc = [6  16  26 ];gal = [91 101];
[GLC,GAL] = meshgrid(glc,gal);
hold on;plot(gal_res(GAL),glc_res(GLC),'kx');
plot(gal_res(GAL),glc_res(GLC),'ko');title(num2str(vmax_hxt))
RR = flipud(R(glc,gal))
FF = flipud(f(glc,gal))

Measured = [0.88 0.86; 0.83 0.78; 0.79 0.7]


%% Draw the c12/c13 uptake experiment
figure(4)

R = (gal_in_hxt'+gal_in_gal2')./(glc_in_hxt'+glc_in_gal2'); % ratio of c12/c13
% R = 1/2^doublings+(2^doublings-1)*R/2^doublings; % one doubling correct for a diffrent growth.
r=R;
f = 1./(1+r)

% r = 4/3*R;


h =surf(log2(agal_vec/55),log2(aglu_vec/55),r);set(h,'edgecolor','none');
view(0,90);xlabel('Gal');ylabel('Glc');title('GAL2 glc');caxis([0 c_max]);

% S glucose -6.5 and galactose -5 ( 3 tubes)
% S glucose -6.5 and galactose -1 ( 3 tubes)    
% S glucose -6.5 and galactose 1 ( 3 tubes)
  
% S glucose -4 and galactose -5 ( 3 tubes)
% S glucose -4 and galactose -1 ( 3 tubes)
% S glucose -4 and galactose 1 ( 3 tubes)
glc = [6  16  26 36];gal = [91 101];
[GLC,GAL] = meshgrid(glc,gal);
hold on;plot(gal_res(GAL),glc_res(GLC),'kx');
plot(gal_res(GAL),glc_res(GLC),'ko');title(num2str(vmax_hxt))
RR = flipud(r(glc,gal))
FF = flipud(f(glc,gal))

