function dy = MurrayEqs(t,y,agal)



%% Parameters
kf81 = 100;
kr81 = 1500;
kf83 = 100;
kr83 = 1;
kf84 = 100;
kr84 = 25;

aG1 = 15;
aG3 = 0.9;
aG4 = 0.2;
a0G80 = 0.6;
aG80 = 0.9;

KG1 = 8;
KG3 = 8;
KG80 = 2;

n1 = 3;
n3 = 2;
n80 = 2;

gG1 = 0.004;
gG3 = 0.004;
gG4 = 0.004;
gG80 = 0.004;
gG81 = 0.004;
gG83 = 0.004;
gG84 = 0.004;
 
epsilon = 0.1;

aG1s = 0.1;
aG3s = 0.1;
aG80s = 1.5;

%% Effctive paramters
w = kr81*kf81/( kr81 + kf81 ) - kf81;
d = kr83*kf83/( kr83 + kf83 ) - kf83;
b = kr84*kf84/( kr84 + kf84 ) - kf84;
%%
% y(1) = G1;
% y(2) = G80;
% y(3) = G3;
% y(4) = G4;
for i = 1:4
    y(i) = max([0 y(i)]);
end
dy = zeros(4,1);    % a column vector
dy(1) = agal*epsilon + aG1*(y(4)^n1/(KG1^n1 + y(4)^n1)) + w*y(1)*y(2) - gG1*y(1);
dy(2) = a0G80 + aG80*(y(4)^n80/(KG80^n80 + y(4)^n80)) + w*y(1)*y(2) + d*y(3)*y(2) + b*y(4)*y(2) - gG80*y(2);
dy(3) = agal + aG3*(y(4)^n3/(KG1^n3 + y(4)^n3)) + d*y(3)*y(2) - gG3*y(3);
dy(4) = aG4 + b*y(4)*y(2) - gG4*y(4);

