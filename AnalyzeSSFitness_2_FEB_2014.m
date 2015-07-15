%% Analyze Facs




clear all
close all
clc

home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\';
cd(dir_name);
load gate
cd(home_name);



home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\t0\';
cd(dir_name);
filename = 'SteadyStateFitness_2_FEB_2014_t0'
load(filename)
cd(home_name);

%%
fsc_ind = 2;
ssc_ind = 5;

yfp_ind = 3;
rfp_ind = 6;
bfp_ind = 8;

well_num_str= {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12' };
well_letter = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

strain_wells = {'A01' 'A03' 'A05' 'A07' 'A09''A011'};
strain_name = {'wo1' 'wo2' 'wo3' 'w1' 'w2' 'w3'};

ctrs = logspace(-4,4,200);

%%
for i = 1:6

            
            M_yfp(i,:) = hist(data{i}(:,yfp_ind),ctrs);
            M_rfp(i,:)  = hist(data{i}(:,rfp_ind),ctrs);
            M_bfp(i,:)  = hist(data{i}(:,bfp_ind),ctrs);

            FSC{i}  = data{i}(:,fsc_ind);
            SSC{i}  =  data{i}(:,ssc_ind);
            
            counts(i) = length(FSC{i});
            

 end

%%
close all

figure(1)
subplot(3,1,1)
plot(ctrs,M_yfp');set(gca,'xscale','log');title('YFP')
figure(1)
subplot(3,1,2)
plot(ctrs,M_rfp');set(gca,'xscale','log');title('rFP')
figure(1)
subplot(3,1,3)
plot(ctrs,M_bfp');set(gca,'xscale','log');title('bFP')

% %% Segment
% close all
% 
% for i = 6:6
%     figure(i)
%     plot((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),'.k','markersize',4);hold on
%     set(gca,'xscale','log','yscale','log');xlabel('rfp');ylabel('bfp');
%     % CROP
%     title('RFP')
%    [x_rfp,y_rfp]=getline(gcf,'closed')
%     title('BFP')
%     [x_bfp,y_bfp]=getline(gcf,'closed')
%     
%     ind_red = find(inpolygon((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),x_rfp,y_rfp));
%     ind_blue = find(inpolygon((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),x_bfp,y_bfp));
% 
% % 
% % %     ind_rfp = find((data{i}(:,rfp_ind)>cut_rfp).*((data{i}(:,bfp_ind)<cut_bfp)));
% % %     ind_bfp = find((data{i}(:,rfp_ind)<cut_rfp).*((data{i}(:,bfp_ind)>cut_bfp)));
% % %     ind_d = find((data{i}(:,rfp_ind)>cut_rfp).*((data{i}(:,bfp_ind)>cut_bfp)));
% % %     ind_back = find((data{i}(:,rfp_ind)<cut_rfp).*((data{i}(:,bfp_ind)<cut_bfp)));
% % %     
% % %     s = ones(size((data{i}(:,rfp_ind))));
% % %     s(ind_bfp)=2;
% % %     s(ind_back)=3;
% % %     s(ind_d)=4;
% % %     
% %     options = statset('MaxIter',400);
% % 
% %     obj= gmdistribution.fit([(data{i}(:,rfp_ind)),(data{i}(:,bfp_ind))],2,'Replicates',3,'Options',options);
% %     [idx,nlogl,P,logpdf,M] = cluster(obj,[(data{i}(:,rfp_ind)),(data{i}(:,bfp_ind))]);
% % %     [temp idx_ind] = max( [mean( data{i}(find(idx_rfp==1),rfp_ind)), mean( data{i}(find(idx_rfp==2),rfp_ind)  )])
% %     
% %     ind_1= find((idx==1));
% %     ind_2 = find((idx==2));
% % 
% % 
% % %     obj_bfp = gmdistribution.fit([(data{i}(:,bfp_ind)),(data{i}(:,bfp_ind))],2,'Replicates',3,'Options',options);
% % %     [idx_bfp,nlogl,P,logpdf,M] = cluster(obj,[(data{i}(:,bfp_ind)),(data{i}(:,ssc_ind))]);
% % %     [temp idx_ind] = max( [mean( data{i}(find(idx_bfp==1),bfp_ind)), mean( data{i}(find(idx_bfp==2),bfp_ind)  )])
% % %     ind_b = ((idx_bfp==idx_ind));
% % % 
% % %     ind_red = find(ind_r.*~ind_b);
% % %     ind_blue = find(ind_b.*~ind_r);
% % 
% % 
% %     cutoff = 0.05;
% %     F1 = mvncdf([(data{i}(ind_1,rfp_ind)),(data{i}(ind_1,bfp_ind))],obj.mu(1,:),obj.Sigma(:,:,1));
% %     ind_F1 = find ( (F1>cutoff).*(F1<(1-cutoff)) );
% %     ind_1 = ind_1(ind_F1);
% % 
% %     F2 = mvncdf([(data{i}(ind_2,rfp_ind)),(data{i}(ind_2,bfp_ind))],obj.mu(2,:),obj.Sigma(:,:,2));
% %     ind_F2 = find ( (F2>cutoff).*(F2<(1-cutoff)) );
% %     ind_2 = ind_2(ind_F2);
% %     
%     figure(i);hold on;
%     plot((data{i}(ind_red,rfp_ind)),(data{i}(ind_red,bfp_ind)),'.r','markersize',4);
%     plot((data{i}(ind_blue,rfp_ind)),(data{i}(ind_blue,bfp_ind)),'.b','markersize',4);
%     
% end
%%
close all;
cut_rfp = 12;
cut_bfp = 20;
for i = 1:6
    
%     ind_rfp = find((data{i}(:,rfp_ind)>cut_rfp).*((data{i}(:,bfp_ind)<cut_bfp)));
%     ind_bfp = find((data{i}(:,rfp_ind)<cut_rfp).*((data{i}(:,bfp_ind)>cut_bfp)));
%     ind_d = find((data{i}(:,rfp_ind)>cut_rfp).*((data{i}(:,bfp_ind)>cut_bfp)));
%     ind_back = find((data{i}(:,rfp_ind)<cut_rfp).*((data{i}(:,bfp_ind)<cut_bfp)));
%     
    ind_rfp = find(inpolygon((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),x_rfp,y_rfp));
    ind_bfp = find(inpolygon((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),x_bfp,y_bfp));

    counts_rfp(i) = length(ind_rfp);
    counts_bfp(i) = length(ind_bfp);

    y_b = hist(data{i}(ind_bfp,yfp_ind),ctrs);
    y_r = hist(data{i}(ind_rfp,yfp_ind),ctrs);

    figure(1)
    subplot(2,3,i)
    plot((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),'.k','markersize',4);hold on;
    plot((data{i}(ind_rfp,rfp_ind)),(data{i}(ind_rfp,bfp_ind)),'.r','markersize',4);
    plot((data{i}(ind_bfp,rfp_ind)),(data{i}(ind_bfp,bfp_ind)),'.b','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
%     plot((data{i}(ind_d,rfp_ind)),(data{i}(ind_d,bfp_ind)),'.m','markersize',4);
%     plot((data{i}(ind_back,rfp_ind)),(data{i}(ind_back,bfp_ind)),'.c','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
% 
    xlabel('RFP');ylabel('BFP');
    
    figure(2);
    subplot(1,2,1);hold on
    plot(ctrs,y_b/max(y_b),'b');set(gca,'xscale','log');title('A10')
    xlabel('YFP');xlim([0.1 10^4]);
    
    subplot(1,2,2);hold on
     plot(ctrs,y_r/max(y_r),'r');set(gca,'xscale','log');title('H3')
    xlabel('YFP');xlim([0.1 10^4]);


end

%%
figure(4)

bar([counts_bfp;counts_rfp]')
od_bfp = counts_bfp/400*1000/(3*10^7)
od_rfp = counts_rfp/400*1000/(3*10^7)
figure(5)
plot(od_bfp,'ob');hold on;
plot(od_rfp,'or');hold on;
%%
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\t0'
cd(dir_name);
save t0;
cd(home_name);

%% Compare time points
clear all
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\t0'
cd(dir_name);
load t0;
counts_H3(1,:) = counts_rfp;
counts_A10(1,:) = counts_bfp;

od_H3(1,:) = od_rfp;
od_A10(1,:) = od_bfp;
t0 = cell2mat(time_vec)

dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\t1'
cd(dir_name);
load t1;
counts_H3(2,:) = counts_rfp;
counts_A10(2,:) = counts_bfp;
od_H3(2,:) = od_rfp;
od_A10(2,:) = od_bfp;
t1 = cell2mat(time_vec)


cd(home_name);

%%

figure(10);hold on
plot(counts_H3(:,1:3),'o-r')
plot(counts_H3(:,4:6),'x-r')
plot(counts_A10(:,1:3),'o-b')
plot(counts_A10(:,4:6),'x-b')
ylabel('counts/400ul')

figure(11);hold on
plot(od_H3(:,1:3),'o-r')
plot(od_H3(:,4:6),'x-r')
plot(od_A10(:,1:3),'o-b')
plot(od_A10(:,4:6),'x-b')
ylabel('OD')



figure(12);
plot([1,1,1],181./(-log2(od_A10(1,1:3)./od_A10(2,1:3))),'ob');hold on;
plot([1,1,1],181./(-log2(od_H3(1,1:3)./od_H3(2,1:3))),'or');hold on;

plot([1,1,1]*2,181./(-log2(od_A10(1,4:6)./od_A10(2,4:6))),'ob');hold on;
plot([1,1,1]*2,181./(-log2(od_H3(1,4:6)./od_H3(2,4:6))),'or');hold on;

xlim([0 3])
ylabel('Doubling time min')
ylim([50 300])

figure(13);
plot([1,1,1],(-log(od_A10(1,1:3)./od_A10(2,1:3)))/3,'ob');hold on;
plot([1,1,1],(-log(od_H3(1,1:3)./od_H3(2,1:3)))/3,'or');hold on;

plot([1,1,1]*2,(-log(od_A10(1,4:6)./od_A10(2,4:6)))/3,'ob');hold on;
plot([1,1,1]*2,(-log(od_H3(1,4:6)./od_H3(2,4:6)))/3,'or');hold on;

xlim([0 3])
ylabel('growth rate 1/hr');
ylim([0 0.7])

figure(15)
A10_gal = mean((-log(od_A10(1,1:3)./od_A10(2,1:3)))/3);
eA10_gal = std((-log(od_A10(1,1:3)./od_A10(2,1:3)))/3);
A10_glc = mean((-log(od_A10(1,4:6)./od_A10(2,4:6)))/3);
eA10_glc = std((-log(od_A10(1,4:6)./od_A10(2,4:6)))/3);


H3_gal = mean((-log(od_H3(1,1:3)./od_H3(2,1:3)))/3);
eH3_gal = std((-log(od_H3(1,1:3)./od_H3(2,1:3)))/3);
H3_glc = mean((-log(od_H3(1,4:6)./od_H3(2,4:6)))/3);
eH3_glc = std((-log(od_H3(1,4:6)./od_H3(2,4:6)))/3);

errorbar (1, A10_gal./A10_glc, sqrt( (1/A10_glc*eA10_gal)^2+(A10_gal./A10_glc^2*eA10_glc)^2 ),'o');hold on
errorbar(2, H3_gal./H3_glc,sqrt( (1/H3_glc*eH3_gal)^2+(H3_gal./H3_glc^2*eH3_glc)^2 ),'o')

ylabel('Fitness deacrse');ylim([0 1.2]);title(['H3 = ',num2str(H3_gal./H3_glc),'+- ',num2str(sqrt( (1/H3_glc*eH3_gal)^2+(H3_gal./H3_glc^2*eH3_glc)^2 ))...
    ,'A10 = ',num2str(A10_gal./A10_glc),'+- ',num2str(sqrt( (1/A10_glc*eA10_gal)^2+(A10_gal./A10_glc^2*eA10_glc)^2 ))]);