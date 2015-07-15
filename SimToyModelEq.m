function dy = SimToyModelEq(t,y,agal,aglu,e)



%% Parameters
d = 0.004; % 90 min 

aG40 = 1;
aG80_0 = 1;
aG3 = 1;

gG4 = d;
gMig = d;
gG80 = d;
gG84 = d;
gG3 = d;
gG83 = d;

kfa = 0.1; kra = 1; % Gives K of 10;
kfr = 0.1; krr = 1;

kr84 = 25;
kf84 = 100;

kr83 = 1;
kf83 = 100;

KG3 = 1;
n3 = 2;

T = 1;

%%

for i = 1:11
    y(i) = max([0 y(i)]);
end
dy = zeros(11,1);    % a column vector

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


%% Gal3

dy(3) = agal - gG3*y(3)...
        - kf83*y(3)*y(2) + kr83*y(6);...
%         + aG3*((y(8))^n3/(KG3^n3 + (y(8))^n3));
    
dy(6) = kf83*y(3)*y(2) - kr83*y(6) - gG83*y(6);


%% Gal80 

    
dy(2) = aG80_0 - gG80*y(2)...
            - kf84*y(4)*y(2) + kr84*y(7)...
            - kf83*y(3)*y(2) + kr83*y(6);
% %         + aG80*(y(8))^n80/(KG80^n80 + (y(8))^n80)...
% %         - kf81*y(1)*y(2) + kr81*y(5)...



dy(7) = kf84*y(4)*y(2) - kr84*y(7) - gG84*y(7); 

%% DNA binding and unbinsing: Gal4 mig1 and DNA

% GAl4
dy(4) = aG40 - gG4*y(4)...
        -(T - y(8) - y(9) - y(10)) * y(4) *kfa + (y(8)+y(10))*kra - y(9)*1/e*kfa*y(4)...
        - kf84*y(4)*y(2) + kr84*y(7);
    
% MIG1    
dy(11) = aglu - gMig*y(11)...
        -(T - y(8) - y(9) - y(10)) * y(11) *kfr + (y(9)+y(10))*krr - y(8)*1/e*kfr*y(11);
     
% E:G4
dy(8) = (T - y(8) - y(9) - y(10)) * y(4) *kfa + y(10)*krr...
        -y(8)*( kra + 1/e*kfr*y(11));
% MIG1:E
dy(9) = (T - y(8) - y(9) - y(10)) * y(11) *kfr + y(10)*kra...
        -y(9)*( krr + 1/e*kfa*y(4));
% MIG1:G2
dy(10) = y(8)*1/e*kfr*y(11) + y(9)*1/e*kfa*y(4)...
         -y(10)*(krr + kra);...

    
  
