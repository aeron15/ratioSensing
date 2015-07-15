%% Simulating growth 
close all
clear all
clc
%%

ra1 = 0.6; %1/hr
ra2 = 0.4;

rh1 = 0.55;
rh2 = 0.35;

tha = 0.1; %
thh = 1;

ti = 3; %hr



d = 0.2; % %/[hrXod]
cost = 0.3;
%%
na0 = 1*0.001; % OD
nh0 = 0*0.001;

countera = 0;
counterh = 0;


c10 = 1;
c20 = 0.25;



Y0 = [na0,nh0,c10,c20,countera,counterh];

[T Y] = ode15s(@(t,y) SimulatingGrowthEq(t,y,ra1,ra2,rh1,rh2,tha,thh,ti,d,cost),[0 30],Y0); % Solve ODE

figure(1)

subplot(2,1,1);hold on;
plot(T,Y(:,1),'k');
plot(T,Y(:,2),'r');
subplot(2,1,2);hold on;
plot(T,Y(:,3),'b')
plot(T,Y(:,4),'c')

figure(2)
plot(T,Y(:,5));


%% single 
C1 = logspace(-3,0.3,10);
C2 = logspace(-3,0.3,5);

na0 = 0.001; % OD
nh0 = 0.001;

countera = 0;
counterh = 0;

timaa = 0;
timah = 0;

for i = 1:length(C1)
    for j = 1:length(C2);

        
        c10 = C1(i);
        c20 = C2(j);
        
        Y0 = [na0,0,c10,c20,countera,counterh];

        [T Y] = ode15s(@(t,y) SimulatingGrowthEq(t,y,ra1,ra2,rh1,rh2,tha,thh,ti,d,cost),[0 30],Y0); % Solve ODE

    yield_single_a(i,j) = Y(end,1);
    end
end
%%

figure(3)
subplot(2,1,1)
surf(log2(C1),log2(C2),yield_single_a');xlabel('Glc');ylabel('Gal')



