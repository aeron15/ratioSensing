function dy = SimToyModelEqEff(t,y,agal,aglu,e,...
                               d,amig,aG4,aG80,aG3,aG1,...
                               K4DNA,KmigDNA,...
                               kr84,kf84,...
                               kr83,kf83,...
                               kgal,amig0 )
                               



%% Parameters
% d = 0.004; % 90 min 

gG4 = d;
gMig = d;
gG80 = d;
gG84 = d;
gG3 = d;
gG83 = d;
gG1 = d;

n80 =2;
n3=3;
n1=3;
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


%% Gal1

dy(1) = aG1*AR - gG1*y(1);

%% Gal3

% dy(8) = aG3*AR - gG3*y(8)...
%         -agal/(2^kgal+gal)*kf3*y(8) + kr3*y(3);
 
dy(3) =  agal/(2^kgal+agal)*aG3*AR - gG3*y(3)...
        - kf83*y(3)*y(2) + kr83*y(6);
    
% dy(3) =  agal*kf3*y(8) - kr3*y(3) - gG3*y(3)...
%         - kf83*y(3)*y(2) + kr83*y(6);
         
dy(6) = kf83*y(3)*y(2) - kr83*y(6) - gG83*y(6);


%% Gal80 

    
dy(2) =  aG80*A - gG80*y(2)...
            - kf84*y(4)*y(2) + kr84*y(7)...
            - kf83*y(3)*y(2) + kr83*y(6);
% %         - kf81*y(1)*y(2) + kr81*y(5)...



dy(7) = kf84*y(4)*y(2) - kr84*y(7) - gG84*y(7); 

%% DNA binding and unbinsing: Gal4 mig1 and DNA

% GAl4
dy(4) = aG4*R- gG4*y(4)...
        - kf84*y(4)*y(2) + kr84*y(7);
    
% MIG1    
dy(11) = amig0*amig + amig*aglu/(aglu+0.2) - gMig*y(11);

    
  
