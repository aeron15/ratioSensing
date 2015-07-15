%% Ratio sensor 
clear all
close all
clc

%
g1 = 0.6; %1/hr
g2 = 0.4; %1/hr
d1 = 0.2; %/OD
d2 = 0.2;
cost = 1.2;
r=1;
ti = 3;
od_0 = 0.001;

   
C1 = [2.^(-11:1:1)];
C2 = [2.^(-11:1:1)];

ts_vec = linspace(0,15,20);

%% Loop over concentrations seqential eating
% 

clear ts_opt_sim


od_0 = 0.001;
for j = 1:length(C1);
    for k =1:length(C2);
        T0(j) = 1/g1*log(1 + C1(j)/(d1*od_0));
        clear F
        for i = 1:length(ts_vec);

            ts = ts_vec(i);
            Y0 = [od_0,C1(j),C2(k),0,0];

            [T Y] = ode15s(@(t,y) SimulatingGrowthSequntiolEq(t,y,g1,g2,d1,d2,ts,ti,cost,r),[0 70],Y0); % Solve ODE

%                 figure(1)
%                 plot(T,Y)

            ind = find(abs((Y(:,2)))<0.0001);T1 = T(ind(1));
            ind = find(abs((Y(:,3)))<0.0001);T2 = T(ind(1));
            Tend(i) = max(T1,T2);
            BMend(i) = Y(end,1);
            F(i) = BMend(i)/Tend(i);

        end
        [temp,ind] = max(F);
        ts_opt_sim(j,k) = ts_vec(ind);
    end
end


figure(10)
pcolor(log2(C2),log2(C1),ts_opt_sim);

ind_mat = ts_opt_sim;
ind_mat(find(ts_opt_sim<=0)) = 1;
ind_mat(find(ts_opt_sim>0))=0;

figure(11)
pcolor(log2(C2),log2(C1),ind_mat);

figure(12)
plot(ts_opt_sim);hold on;plot(T0-ti/cost,'k','linewidth',2)



%% Pareael eating
clear ts_opt_sim

for j = 1:length(C1);
    for k =1:length(C2);
        T0(j) = 1/g1*log(1 + C1(j)/(d1*od_0));
        clear F
        for i = 1:length(ts_vec);

            ts = ts_vec(i);
            Y0 = [od_0,C1(j),C2(k),0,0];

            [T Y] = ode15s(@(t,y) SimulatingGrowthParalelEq(t,y,g1,g2,d1,d2,ts,ti,cost,r),[0 70],Y0); % Solve ODE

%                 figure(1)
%                 plot(T,Y)

            ind = find(abs((Y(:,2)))<0.0001);T1 = T(ind(1));
            ind = find(abs((Y(:,3)))<0.0001);T2 = T(ind(1));
            Tend(i) = max(T1,T2);
            BMend(i) = Y(end,1);
            F(i) = BMend(i)/Tend(i);

        end
        [temp,ind] = max(F);
        ts_opt_sim(j,k) = ts_vec(ind);
    end
end


figure(13)
pcolor(log2(C2),log2(C1),ts_opt_sim);

ind_mat = ts_opt_sim;
ind_mat(find(ts_opt_sim<=0)) = 1;
ind_mat(find(ts_opt_sim>0))=0;

figure(14)
pcolor(log2(C2),log2(C1),ind_mat);


figure(15)
plot(ts_opt_sim);hold on;plot(T0-ti/cost,'k','linewidth',2)