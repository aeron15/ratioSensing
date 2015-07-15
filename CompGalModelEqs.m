function dy = CompGalModelEqs(t,y,agal,aglu,e,x)


%% Parameters

kf83 = x(1);
kr83 = x(2);
kf84 = x(3);
kr84 = x(4);
aG1 =  x(5) ;
aG3 =  x(6);
aG4 = x(7);
a0G80 = x(8);
aG80 = x(9);
amig = x(10);
KG1 = x(11);
KG3 = x(12);
KG80 = x(13);
Kmig1 = x(14);
Kmig3 = x(15);
Kmig4 = x(16);

n1 = x(17);
n3 = x(18);
n80 = x(19);

gG1 = x(20);
gG3 = x(21);
gG4 = x(22);
gG80 = x(23);
gG33 = x(24);
gG83 = x(25);
gG84 = x(26);
gMig = x(27);
Kgal = x(28);
Kglu = x(29);

kf33 = x(30);
kr33 = x(31);

% vmax_hxt = x(32);
% Khxt_glc = x(33);
% Khxt_gal = x(34);
% 
% vmax_gal2 = x(35);
% Kgal2_glc = x(36);
% Kgal2_gal = x(37);


%%


for i = 1:11
    y(i) = max([0 y(i)]);
end
dy = zeros(11,1);    % a column vector

%% Gal genes


% Gal 1
dy(1) = aG1* (1/ ( (1 + y(11)/Kmig1)*(1 + (KG1/y(4))^n1) +(1/e-1)*(y(11)/Kmig1)) )  - gG1*y(1);

%G80    
dy(2) = a0G80 - gG80*y(2)...
        + aG80*(y(4))^n80/(KG80^n80 + (y(4))^n80)...
        - kf83*y(5 )*y(2) + kr83*y(6)...
        - kf84*y(4)*y(2) + kr84*y(7);

    
% Gal3
dy(3) =  aG3 * (1/ ( (1 + y(11)/Kmig3)*(1 + (KG3/y(4))^n3) +(1/e-1)*(y(11)/Kmig3)   )) - gG3*y(3)...
         - kf33*agal*y(3) + kr33*y(5);


% Gal3 activated    
dy(5) =    kf33*agal*y(3) - kr33*y(5)...
           - kf83*y(5)*y(2) + kr83*y(6)...
           - gG33*y(5); 
    
% G380    
dy(6) = kf83*y(5)*y(2) - kr83*y(6) - gG83*y(6);
   
%G804
dy(7) = kf84*y(4)*y(2) - kr84*y(7) - gG84*y(7); 



%% DNA binding and unbinsing: Gal4 mig1 and DNA
% GAl4
dy(4) = aG4*1/(1+y(11)/Kmig4) - gG4*y(4)...
        - kf84*y(4)*y(2) + kr84*y(7);
    
% MIG1    
dy(11) = amig*aglu^3/(aglu+Kglu)^3 - gMig*y(11);
     
  
