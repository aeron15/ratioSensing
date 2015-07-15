%% Analyze Facs

clear all
close all
clc


home_name = pwd;
dir_name = 'C:\Users\ys151\Desktop\WorkingZone\Gal\GlcAnalogFacs21MAR2014\HighIntensity\';
cd(dir_name);
filename = 'GlcAnalogFacs21MAR2014HighIntensity'
load(filename)
cd(home_name);

%%
fsc_ind = 1;
ssc_ind = 2;

yfp_ind = 3;
rfp_ind = 6;



ctrs = logspace(-4,4,200);
%%

well_num_str= {'01' '02' '03' '04' '05' '06'};
well_letter = {'A' 'B'};

strain_wells = {'A01' 'A03' 'A05' 'A07' 'B02' 'B04' 'C01' 'C03' 'C05'};
strain_name = {'M1' 'M2' 'M3' 'M4' 'M5' 'M6' 'M7' 'M8' };

ctrs = logspace(-4,4,200);

%%
for i = 1:length(well_num_str)
    k=1;
    i
    for j = 1:length(data)
        x=cell2mat(well_id_plate{j});
        if find(strcmp(x(2:3),well_num_str))==i
           
            
            M_rfp{i,k} = hist(data{j}(:,rfp_ind),ctrs);
            M_gfp{i,k} = hist(data{j}(:,yfp_ind),ctrs);
            
            x_rfp(i,k,:) = (data{j}(1:50000,rfp_ind));
            x_gfp(i,k,:) = (data{j}(1:50000,yfp_ind));

            FSC{i,k} = data{j}(:,fsc_ind);
            SSC{i,k} =  data{j}(:,ssc_ind);
            
            counts{i,k} = length(FSC{i,k});
            
            time{i,k} = time_vec{j};
           k=k+1 
        end
    end
end
%%
figure(1)
 subplot(2,2,1)
 pcolor([1:6],ctrs,cell2mat(M_gfp(:,1))');set(gca,'yscale','linear')
  subplot(2,2,2)
   pcolor([1:6],ctrs,cell2mat(M_rfp(:,1))');set(gca,'yscale','linear')

 subplot(2,2,3)
  pcolor([1:6],ctrs,cell2mat(M_gfp(:,2))');set(gca,'yscale','linear')

 subplot(2,2,4)
 pcolor([1:6],ctrs,cell2mat(M_rfp(:,2))');set(gca,'yscale','linear')

%%
figure(2)
 subplot(2,2,1)
h = boxplot((squeeze(x_gfp(:,1,:)))');ylim([0 150]);set(h(7,:),'Visible','off');

subplot(2,2,2)
h = boxplot((squeeze(x_gfp(:,2,:)))');ylim([0 150]);
set(h(7,:),'Visible','off');

 subplot(2,2,3)
h = boxplot((squeeze(x_rfp(:,1,:)))');set(h(7,:),'Visible','off');
ylim([0 350]);

 subplot(2,2,4)
h = boxplot((squeeze(x_rfp(:,2,:)))');set(h(7,:),'Visible','off');
ylim([0 350]);

 