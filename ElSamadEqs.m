function dy = ElSamadEqs(t,y,gal,glu,e)

%% Parameters


kf_gal_g1 = 70.3;      
kr_gal_g1 = 3391;      
kf_glc_mig1 = 56.2;    
kr_glc_mig1 = 3564;

kf81 = 41.1;
kr81 = 700.1;
kf84 = 95.2;
kr84 = 1237;

a0g1 = 0.93;
ag1 = 16.5;
ag4 = 2;
a0g80 = 0.71;
ag80 = 2;
amig1 = 0.75;

Kg1 = 41.6;
K_mig1_g1 = 67.4;
K_mig1_g4 = 33.8;
Kg80 = 14;

ng1 = 3;
n_mig1_g1 = 2;
n_mig1_g4 = 1;
n80 = 2;

dg1 = 0.005;
dg1s = 0.087;
dg4 = 0.004;
dg80 = 0.005;
dr = 0.009;
drs = 0.009;

d81 = 0.015;
d84 = 0.024;



for i = 1:8
    y(i) = max([0 y(i)]);
end
dy = zeros(8,1);    % a column vector


% G1:
dy(1) = a0g1 + ag1*(y(4)^ng1/(y(4)^ng1 + Kg1^ng1))*(K_mig1_g1^n_mig1_g1/(K_mig1_g1^n_mig1_g1 + y(8)^n_mig1_g1)) - kf_gal_g1*gal*y(1) + kr_gal_g1*y(2) -y(1)*dg1;

% G1*:
dy(2) = kf_gal_g1*gal*y(1) - kr_gal_g1*y(2) -kf81*y(2)*y(3) + kr81*y(5) -y(2)*dg1s;
    
%G801:    
dy(5) = kf81*y(2)*y(3) - kr81*y(5) - y(5)*d81;
    
% G4:

dy(4) = ag4*(K_mig1_g4^n_mig1_g4/(K_mig1_g4^n_mig1_g4 + y(8)^n_mig1_g4)) - kf84*y(4)*y(3) + kr84*y(6) - y(4)*dg4;
 
    
%G804:    
dy(6) = kf84*y(4)*y(3)  - kr84*y(6) - y(6)*d84;


%G80
dy(3) = a0g80 + ag80*(y(4)^n80/(y(4)^n80 + Kg80^n80)) - kf81*y(2)*y(3) + kr81*y(5) -kf84*y(4)*y(3) + kr84*y(6) -y(3)*dg80;
      
% mig1
dy(7) = amig1 - kf_glc_mig1*y(7)*glu + kr_glc_mig1*y(8) -y(7)*dr;

% mig 1 activated
dy(8) = kf_glc_mig1*y(7)*glu - kr_glc_mig1*y(8) -y(8)*drs;


 