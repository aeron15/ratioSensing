%% Analyze Facs

clear all
close all
clc


home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\LowOdGrowth_7_APR_2014\';
cd(dir_name);
filename = 'LowOdGrowth_7_APR_2014'
load(filename)
cd(home_name);

%%
fsc_ind = 2;
ssc_ind = 5;

yfp_ind = 3;
rfp_ind = 6;
bfp_ind = 8;
time_ind = 10;


ctrs = logspace(-4,4,200);
%%

well_num_str= {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well_letter = {'A' 'B' 'C' 'D'};
% 
% strain_wells = {'A01' 'A03' 'A05' 'A07' 'B02' 'B04' 'C01' 'C03' 'C05'};
% strain_name = {'M1' 'M2' 'M3' 'M4' 'M5' 'M6' 'M7' 'M8' };

ctrs = logspace(-4,4,200);

%%
for i = 1:length(well_letter)
    k=1;
    i
    for j = 1:length(data)
        x=cell2mat(well_id_plate{j});
        if find(strcmp(x(1),well_letter))==i
           
            
            M_rfp{i,k} = hist(data{j}(:,rfp_ind),ctrs);
            
            FSC{i,k} = data{j}(:,fsc_ind);
            SSC{i,k} =  data{j}(:,ssc_ind);
            
            counts{i,k} = length(FSC{i,k});
            well{i,k} =  well_id_plate{j};
            med_gfp{i,k} = median(data{j}(:,rfp_ind));
            time{i,k} = time_vec{j};
           k=k+1 
        end
    end
end
%%
time = time(:,1:12);
 d_min = datevec(min(min(cell2mat(time))));
 [Nr,Nt] = size(time);
 
 for i = 1:Nr
     
     d = datevec([time{i,:}]);
     t(i,:) = etime(d,repmat(d_min,[length(d) 1]))/3600;
 end
 %%
% for i = 1:length(well_num_str)
%     
%     M =  cell2mat({M_rfp{i,:}}');
%     d = datevec([time{i,:}]);
%     t_strain = [0;cumsum(etime(d(2:end,:),d(1:end-1,:)))/3600];
% 
%     figure(1)
%     subplot(3,3,i)
%     surf(t_strain,ctrs,M','FaceColor','flat','EdgeColor','none');
%     set(gca,'yscale','linear');grid off;box off;axis tight;  view(0,90);
%     title(strain_name{i});
%     xlim([0 6]);ylim([0 1000])
%      
% end
% Set_fig_YS(figure(1),12,12,12)
%%
counts = counts(:,1:12);
counts_mat = cell2mat(counts(:,1:12));
figure(2)
plot(t,counts_mat,'.-');
set(gca,'yscale','linear');

%%

for i = 1:12
    
    s = fit(t(:,i),log(counts_mat(:,i)),'a*x+b','robust','off');
    m(i) = s.a;
    C(:,:,i) = confint(s,0.682);
    figure(5)
    subplot(4,3,i);
    plot(t(:,i),log(counts_mat(:,i)),'.');hold on;
    plot(t(:,i),s(t(:,i)));
    [R,P] = corrcoef(s(t(:,i)),log(log(counts_mat(:,i))));
    title(['m = ',num2str(s.a),'\rho = ',num2str(R(1,2)),' p = ',num2str(P(1,2))]);
end
figure(4)
eu = squeeze(C(2,1,:))';
ed = squeeze(C(1,1,:))';

errorbar([1:12],m,m-ed,eu-m,'o');

