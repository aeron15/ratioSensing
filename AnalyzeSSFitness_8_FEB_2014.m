%% Analyze Facs

%%%%%%%%%%%%%%%%%%%%%%%%%% WT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\';
cd(dir_name);
load gate
cd(home_name);


home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_8_FEB_2014\wt\';
cd(dir_name);
filename = 'SteadyStateFitness_8_FEB_2014_wt'
load(filename)
cd(home_name);

%%
fsc_ind = 2;
ssc_ind = 5;

yfp_ind = 3;
rfp_ind = 6;
bfp_ind = 8;



ctrs = logspace(-4,4,200);

%% Segment
close all;
cut_rfp = 12;
cut_bfp = 20;
for i = 1:18
    
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
    subplot(5,5,i)
    plot((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),'.k','markersize',4);hold on;
    plot((data{i}(ind_rfp,rfp_ind)),(data{i}(ind_rfp,bfp_ind)),'.r','markersize',4);
    plot((data{i}(ind_bfp,rfp_ind)),(data{i}(ind_bfp,bfp_ind)),'.b','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
%     plot((data{i}(ind_d,rfp_ind)),(data{i}(ind_d,bfp_ind)),'.m','markersize',4);
%     plot((data{i}(ind_back,rfp_ind)),(data{i}(ind_back,bfp_ind)),'.c','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
% 
    xlabel('RFP');ylabel('BFP');
    
    figure(2);
    subplot(5,5,i);hold on
    plot(ctrs,y_b/max(y_b),'b');set(gca,'xscale','log');
    
    plot(ctrs,y_r/max(y_r),'r');set(gca,'xscale','log');
    xlabel('YFP');xlim([0.1 10^4]);title(num2str(i));


end

%% divide in conctrations

ind_1_0 = [1:3];
ind_2_0 = [16:18];

ind_1_065 = [7:9];
ind_2_065 = [13:15];

ind_1_2 = [4:6];
ind_2_2 = [10:12];


gal0_wt(1,:) = counts_bfp(ind_1_0);
gal0_wt(2,:) = counts_bfp(ind_2_0);
gal0_H3(1,:) = counts_rfp(ind_1_0);
gal0_H3(2,:) = counts_rfp(ind_2_0);
t(1,:) =  etime(datevec(cell2mat(time_vec(ind_2_0))),datevec(cell2mat(time_vec(ind_1_0))))/3600;

gal065_wt(1,:) = counts_bfp(ind_1_065);
gal065_wt(2,:) = counts_bfp(ind_2_065);
gal065_H3(1,:) = counts_rfp(ind_1_065);
gal065_H3(2,:) = counts_rfp(ind_2_065);
t(2,:) =  etime(datevec(cell2mat(time_vec(ind_2_065))),datevec(cell2mat(time_vec(ind_1_065))))/3600;

gal2_wt(1,:) = counts_bfp(ind_1_2);
gal2_wt(2,:) = counts_bfp(ind_2_2);
gal2_H3(1,:) = counts_rfp(ind_1_2);
gal2_H3(2,:) = counts_rfp(ind_2_2);
t(3,:) =  etime(datevec(cell2mat(time_vec(ind_2_2))),datevec(cell2mat(time_vec(ind_1_2))))/3600;


%% PLOT
close  all

figure(1)

figure(1);hold on;
plot([1,1,1],(-log(gal0_wt(1,:)./gal0_wt(2,:)))/t(1,1),'.b');hold on;
plot([1,1,1]*2,(-log(gal065_wt(1,:)./gal065_wt(2,:)))/t(2,1),'.b');hold on;
plot([1,1,1]*3,(-log(gal2_wt(1,:)./gal2_wt(2,:)))/t(3,1),'.b');hold on;

plot([1,1,1],(-log(gal0_H3(1,:)./gal0_H3(2,:)))/t(1,1),'.r');hold on;
plot([1,1,1]*2,(-log(gal065_H3(1,:)./gal065_H3(2,:)))/t(2,1),'.r');hold on;
plot([1,1,1]*3,(-log(gal2_H3(1,:)./gal2_H3(2,:)))/t(3,1),'.r');hold on;

plot([1 2 3],[median((-log(gal0_wt(1,:)./gal0_wt(2,:)))/t(1,1)),...
              median((-log(gal065_wt(1,:)./gal065_wt(2,:)))/t(2,1)),...
              median((-log(gal2_wt(1,:)./gal2_wt(2,:)))/t(3,1))],'-ob')

          plot([1 2 3],[median((-log(gal0_H3(1,:)./gal0_H3(2,:)))/t(1,1)),...
              median((-log(gal065_H3(1,:)./gal065_H3(2,:)))/t(2,1)),...
              median((-log(gal2_H3(1,:)./gal2_H3(2,:)))/t(3,1))],'-or')

%%
figure(2);hold on;
plot(gal0_wt,'-ro')
plot(gal065_wt,'-bo')
plot(gal2_wt,'-ko')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%A10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A10

clear all
close all
clc

home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_2_FEB_2014\';
cd(dir_name);
load gate
cd(home_name);


home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\SteadyStateFitness_8_FEB_2014\a10\';
cd(dir_name);
filename = 'SteadyStateFitness_8_FEB_2014_a10'
load(filename)
cd(home_name);

%%
fsc_ind = 2;
ssc_ind = 5;

yfp_ind = 3;
rfp_ind = 6;
bfp_ind = 8;



ctrs = logspace(-4,4,200);

%% Segment
close all;
cut_rfp = 12;
cut_bfp = 20;
for i = 1:18
    
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
    subplot(5,5,i)
    plot((data{i}(:,rfp_ind)),(data{i}(:,bfp_ind)),'.k','markersize',4);hold on;
    plot((data{i}(ind_rfp,rfp_ind)),(data{i}(ind_rfp,bfp_ind)),'.r','markersize',4);
    plot((data{i}(ind_bfp,rfp_ind)),(data{i}(ind_bfp,bfp_ind)),'.b','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
%     plot((data{i}(ind_d,rfp_ind)),(data{i}(ind_d,bfp_ind)),'.m','markersize',4);
%     plot((data{i}(ind_back,rfp_ind)),(data{i}(ind_back,bfp_ind)),'.c','markersize',4);set(gca,'xscale','log','yscale','log');hold on;
% 
    xlabel('RFP');ylabel('BFP');
    
    figure(2);
    subplot(5,5,i);hold on
    plot(ctrs,y_b/max(y_b),'b');set(gca,'xscale','log');
    
    plot(ctrs,y_r/max(y_r),'r');set(gca,'xscale','log');
    xlabel('YFP');xlim([0.1 10^4]);title(num2str(i));


end

%% divide in conctrations

ind_1_0 = [6:8];
ind_2_0 = [16:18];

ind_1_065 = [4:5];
ind_2_065 = [13:14];

ind_1_2 = [1:3];
ind_2_2 = [10:12];


gal0_wt(1,:) = counts_bfp(ind_1_0);
gal0_wt(2,:) = counts_bfp(ind_2_0);
gal0_H3(1,:) = counts_rfp(ind_1_0);
gal0_H3(2,:) = counts_rfp(ind_2_0);
t(1,:) =  etime(datevec(cell2mat(time_vec(ind_2_0))),datevec(cell2mat(time_vec(ind_1_0))))/3600;

gal065_wt(1,:) = counts_bfp(ind_1_065);
gal065_wt(2,:) = counts_bfp(ind_2_065);
gal065_H3(1,:) = counts_rfp(ind_1_065);
gal065_H3(2,:) = counts_rfp(ind_2_065);
t(2,1:2) =  etime(datevec(cell2mat(time_vec(ind_2_065))),datevec(cell2mat(time_vec(ind_1_065))))/3600;

gal2_wt(1,:) = counts_bfp(ind_1_2);
gal2_wt(2,:) = counts_bfp(ind_2_2);
gal2_H3(1,:) = counts_rfp(ind_1_2);
gal2_H3(2,:) = counts_rfp(ind_2_2);
t(3,:) =  etime(datevec(cell2mat(time_vec(ind_2_2))),datevec(cell2mat(time_vec(ind_1_2))))/3600;

%%
%% PLOT
close  all

figure(1)

figure(1);hold on;
plot([1,1,1],(-log(gal0_wt(1,:)./gal0_wt(2,:)))/t(1,1),'.b');hold on;
plot([1,1,1]*2,(-log(gal065_wt(1,:)./gal065_wt(2,:)))/t(2,1),'.b');hold on;
plot([1,1,1]*3,(-log(gal2_wt(1,:)./gal2_wt(2,:)))/t(3,1),'.b');hold on;

plot([1,1,1],(-log(gal0_H3(1,:)./gal0_H3(2,:)))/t(1,1),'.r');hold on;
plot([1,1,1]*2,(-log(gal065_H3(1,:)./gal065_H3(2,:)))/t(2,1),'.r');hold on;
plot([1,1,1]*3,(-log(gal2_H3(1,:)./gal2_H3(2,:)))/t(3,1),'.r');hold on;

plot([1 2 3],[median((-log(gal0_wt(1,:)./gal0_wt(2,:)))/t(1,1)),...
              median((-log(gal065_wt(1,:)./gal065_wt(2,:)))/t(2,1)),...
              median((-log(gal2_wt(1,:)./gal2_wt(2,:)))/t(3,1))],'-ob')

          plot([1 2 3],[median((-log(gal0_H3(1,:)./gal0_H3(2,:)))/t(1,1)),...
              median((-log(gal065_H3(1,:)./gal065_H3(2,:)))/t(2,1)),...
              median((-log(gal2_H3(1,:)./gal2_H3(2,:)))/t(3,1))],'-or')

%%
figure(2);hold on;
plot(gal0_wt,'-ro')
plot(gal065_wt,'-bo')
plot(gal2_wt,'-ko')


