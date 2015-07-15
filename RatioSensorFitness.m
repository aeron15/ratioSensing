%% Ratio sensor 
clear all
close all
clc

%
g1 = 0.6; %1/hr
g2 = 0.4; %1/hr
d1 = 0.2; %/OD
d2 = 0.2;
ti = 2;
od_0 = 0.001;
C1 = [0 2.^(-7:1)];
C2 = [0 2.^(-9:1)];
%

cost = 1.2;%linspace(1.1,1.5,5);
r = 1;
ti  = linspace(1,2,3);
ts_vec = linspace(0,20,100);

C1 = 0.1; %[0 2.^(-11:-9)];
C2 = 0.1; %[0 2.^(-1:1)];
% use the analytical expressions
% for j = 1:length(ti);
%     for k =1:length(cost);
%         T0 = 1/g1*log(1 + C1/(d1*od_0));
%         ts_vec = linspace(0,T0,200);
%         for i = 1:length(ts_vec);
%             
%             
%             cts = C1 - od_0*d1*(exp(g1*ts_vec(i))-1);
%             nts = od_0*exp(g1*ts_vec(i));
%             T1  = cost(k)/g1*log(1 + cts/(d1*r*nts) );
%             nT1 = nts*exp(g1/cost(k)*T1);
%             T2  = 1/g2*log(1 + C2/(d2*nT1));
%             BM(i) = nT1*exp(g2*T2);
%             
%             if ti(j)<=T1
%                 T(i) = ts_vec(i)+T1+T2;
%             else
%                 T(i) = ts_vec(i)+ti(j)+T2;
%             end
%         end
%         [temp,ind] = max(log(BM/od_0)./T);
%         ts_opt(j,k) = ts_vec(ind);
%         T_opt(j,k) = T(ind);
%         BM_opt(j,k) = BM(ind);
%     end
% end
% figure(1)
% plot(ts_vec,T,'-o');hold on;
% plot(ts_vec,BM,'-or')
% plot(ts_vec,20*BM./T,'-og')
%%
% Use simulation


for j = 1:length(ti)
    for k =1:length(cost)
            Y0 = [od_0,C1,C2,0,0];
                        
            T0 = 1/g1*log(1 + C1/(d1*od_0));
            
            ts_vec = linspace(0,T0,200);
            for i = 1:length(ts_vec);

            ts = ts_vec(i);


            [T Y] = ode15s(@(t,y) SimulatingGrowthSequntiolEq(t,y,g1,g2,d1,d2,ts,ti(j),cost(k),r),[0 50],Y0); % Solve ODE
%             %
%                 figure(1)
%                 plot(T,Y)

            ind = find(abs((Y(:,2)))<0.0001);T1 = T(ind(1));
            ind = find(abs((Y(:,3)))<0.0001);T2 = T(ind(1));
            Tend(i) = max(T1,T2);
            BMend(i) = Y(end,1);
            F(i) = BMend(i)/Tend(i);

        end
        [temp,ind] = max(F);
        ts_opt_sim(j,k) = ts_vec(ind)

    end
    [j,k]
end

  %%
T0 = 1/g1*log(1+C1/(d1*od_0));
figure(5)
[x,y] = meshgrid(ti,cost);
ts_opt_theory = T0 - x./y;
% plot(ti,ts_opt','.-');hold all;
plot(ti,ts_opt_theory,'k');hold on
plot(ti,ts_opt_sim,'or');plot(ti,ts_opt_sim,'.r');
refline(0,T0)      
         
%% Loop over concentrations
% 

clear ts_opt_sim
ti = 2;

C1 = [2.^(-7:0.5:0)];
C2 = [2.^(-9:0.5:1)];
od_0 = 0.001;
for j = 1:length(C1);
    for k =1:length(C2);
        T0 = 1/g1*log(1 + C1(j)/(d1*od_0));
        ts_vec = linspace(0,T0,200);
        clear F
        for i = 1:length(ts_vec);

            ts = ts_vec(i);
            Y0 = [od_0,C1(j),C2(k),0,0];

            [T Y] = ode15s(@(t,y) SimulatingGrowthSequntiolEq(t,y,g1,g2,d1,d2,ts,ti,cost,r),[0 50],Y0); % Solve ODE
%             %
%                 figure(1)
%                 plot(T,Y)

            ind = find(abs((Y(:,2)))<0.0001);T1 = T(ind(1));
            ind = find(abs((Y(:,3)))<0.0001);T2 = T(ind(1));
            Tend(i) = max(T1,T2);
            BMend(i) = Y(end,1);
            F(i) = BMend(i)/Tend(i);

        end
        [temp,ind] = max(F);
        ts_opt_sim(j,k) = ts_vec(ind)
    end
end


figure(10)
pcolor(log2(C2),log2(C1),ts_opt_sim);

ind_mat = ts_opt_sim;
ind_mat(find(ts_opt_sim<=0)) = 1;
ind_mat(find(ts_opt_sim>0))=0;

figure(11)
pcolor(log2(C2),log2(C1),ind_mat);




%%


% 
% for i = 1:length(C2)
%     
%     for j = 1:length(C1)
%         
%         
%         Cts = C1 - d*od_0*(exp(g1*ts(j))-1);
%              
%         od_ts = od_0*exp(g1*ts(i));
%         
%         T1 = 1/(g1*c)*log( Cts*log(2)/(d*r*od_ts)+1    );
%         
%         od_T1(i,j) = od_ts*exp(g1*c*(T1));
%         
%         T2(i,j) = 1/(g2)*log(    C2(i)*log(2)/(d*od_T1(i,j))+1      );
%         
%         if ti<=T1
%             T(i,j) = ts(i) + T1 + T2(i,j);
%         else
%             T(i,j) = ts(i) + ti + T2(i,j);
%         end
%         
%         F(i,j) = BM(i,j)./T(i,j);
%     end
% end
% 
% figure(1)
% surf(ts,log2(C2),F)
% 
% figure(2)
% surf(ts,log2(C2),BM)
% 
% figure(3)
% surf(ts,log2(C2),T);view(0,90)
% 
% [val,tt] = min(T')
% % 
% % 
% % 
% % t'
% % %%
% % close all
% % 
% % C = [0:0.05:1];
% % ti = 10;%[0:0.5:8];
% % 
% % R1 = [0 2.^(-7:1)];
% % R2 = [0 2.^(-9:1)];
% % OD = 0.01;%[logspace(-3,1,15)];
% % y1 = R1/d*log(2);
% % y2 = R2/d*log(2);
% % T0 = 1/g1*log(1+y1/OD(1));
% % 
% % for i = 1:length(T0)
% %     for j = 1:length(C);
% %         ind(i,j) = T0(i) - C(j)*ti;
% %     end
% % end
% % 
% % 
% % surf(T0,C,ind');hold on;plot(T0,T0./ti)
% % %%
% % 
% % 
% % 
% % for i = 1:length(R1)
% %     for j = 1:length(R2)
% %         
% %         
% %         T1(i,j) = 1/g1*log(y1(i)/od_0+1);
% %         T2(i,j) = 1/g2*log(y2(j)/(od_0+y1(i))+1);
% %         
% %         F1(i,j) = (y1(i)+y2(j))/(T1(i,j)+T2(i,j)+ti);
% %         
% %         if ti>T1(i,j)/c
% %             F2(i,j) = (y1(i)+y2(j))/(T2(i,j)+ti);
% %         else
% %             F2(i,j) = (y1(i)*ti*c/T1(i,j) + y2(j))/(T2(i,j)+ti);
% %         end
% %         
% %         if F2(i,j)>F1(i,j)
% %             ind(i,j) = 1;
% %         else
% %             ind(i,j) = 0;
% %         end
% %     end
% % end
% % 
% % pcolor(log2(R2),log2(R1),ind)
% %         
% 
