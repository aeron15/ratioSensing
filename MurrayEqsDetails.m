function dy = MurrayEqsDetails(t,y,agal,aglu,e)



%% Parameters
kf81 = 100*0;
kr81 = 1500*0;
kf83 = 100;
kr83 = 1;
kf84 = 100;
kr84 = 500;

aG1 = 15;
aG3 = 10;
aG4 = 0.2;
a0G80 = 0.5*1;
aG80 = 1*1;
amig = 1;

KG1 = 8;
KG3 = 8;
KG80 = 2;

Kmig1 = 10;
Kmig3 = 5;
Kmig4 = 70;

n1 = 3;
n3 = 2;
n80 = 2;

d = 0.004;
gG1 = d;
gG3 = d;
gG4 = d;
gG80 = d;
gG81 = d;
gG83 = d;
gG84 = d;
gMig = 0.004;

epsilon_G3 = 100;
epsilon_G1 = 0.1;

aG1s = 0.1;
aG3s = 0.1;
aG80s = 1.5;

Kgal = (2^-4);
Kglu = (2^-1);
%% Effctive paramters
w = kr81*kf81/( kr81 + kf81 ) - kf81;
d = kr83*kf83/( kr83 + kf83 ) - kf83;
b = kr84*kf84/( kr84 + kf84 ) - kf84;
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

for i = 1:11
    y(i) = max([0 y(i)]);
end
dy = zeros(11,1);    % a column vector


% Gal3 and GAl1
%G1:
dy(1) = 0*agal*epsilon_G1 ...
        + aG1*1/ ( (1 + y(11)/Kmig1)*(1 + (KG1/y(4))^n1) +(1/e-1)*(y(11)/Kmig1)   )...
        - kf81*y(1)*y(2) + kr81*y(5) - gG1*y(1);
%G801:    
dy(5) = kf81*y(1)*y(2) - kr81*y(5) - gG81*y(5);
%G3
dy(3) = - kf83*y(3)*y(2) + kr83*y(6)- gG3*y(3) + agal/(agal+Kgal) * aG3 *10* (1/ ( (1 + y(11)/Kmig3)*(1 + (KG3/y(4))^n3) +(1/e-1)*(y(11)/Kmig3)   )) ;
% G380    
dy(6) = kf83*y(3)*y(2) - kr83*y(6) - gG83*y(6);
   
%G80    
dy(2) = a0G80 - gG80*y(2)...
        + aG80*(y(4))^n80/(KG80^n80 + (y(4))^n80)...
        - kf81*y(1)*y(2) + kr81*y(5)...
        - kf83*y(3)*y(2) + kr83*y(6)...
        - kf84*y(4)*y(2) + kr84*y(7);
%G804
dy(7) = kf84*y(4)*y(2) - kr84*y(7) - gG84*y(7); 
     

%% DNA binding and unbinsing: Gal4 mig1 and DNA

% GAl4
dy(4) = aG4*1/(1+y(11)/Kmig4) - gG4*y(4)...
        - kf84*y(4)*y(2) + kr84*y(7);
    
% MIG1    
dy(11) = amig*aglu/(aglu+Kglu) - gMig*y(11);
     
% % E:G4
% dy(8) = (T - y(8) - y(9) - y(10)) * y(4) *kfa + y(10)*krr...
%         -y(8)*( kra + 1/e*kfr*y(11));
% % MIG1:E
% dy(9) = (T - y(8) - y(9) - y(10)) * y(11) *kfr + y(10)*kra...
%         -y(9)*( krr + 1/e*kfa*y(4));
% % MIG1:G2
% dy(10) = y(8)*1/e*kfr*y(11) + y(9)*1/e*kfa*y(4)...
%          -y(10)*(krr + kra);...

    
  
