%% Analyze Facs

clear all
close all
clc

%%
fsc_ind = 1;
ssc_ind = 2;

yfp_ind = 3;
rfp_ind = 4;
time_ind = 6;

ctrs = logspace(-4,4,200);

well_num_str= {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
well_letter = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
% 
% strain_wells = {'A01' 'A03' 'A05' 'A07' 'B02' 'B04' 'C01' 'C03' 'C05'};
% strain_name = {'M1' 'M2' 'M3' 'M4' 'M5' 'M6' 'M7' 'M8' };

time_int = [50 150];


%% Loop over time points
dir_name = {'\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Stratedigm\GEV_10Jun\6hrs\'}
          
filename = {'6hrs'};

for r = 1:1
    clear med_gfp counts
    load([dir_name{r},filename{r}])
    
    for i = 1:length(well_letter)
        k=1;
        i
        for j = 1:length(data)
            x=cell2mat(well_id_plate{j});
            if find(strcmp(x(1),well_letter))==i
                
                
                M_rfp{i,k} = hist(data{j}(:,rfp_ind),ctrs);
                
                FSC{i,k} = data{j}(:,fsc_ind);
                SSC{i,k} =  data{j}(:,ssc_ind);
                rate{i,k} = data{j}(:,time_ind);
                
                x = rate{i,k};rate_avg(i,k) = length((x(x>time_int(1)&x<time_int(2))))/diff(time_int);
                
                counts{i,k} = length(FSC{i,k});
                well{i,k} =  well_id_plate{j};
                med_gfp{i,k} = median(data{j}(:,yfp_ind));
                time{i,k} = time_vec{j};
                
                k=k+1
            end
        end
    end
    
    counts_time{r} = (counts);
    med_gfp_time{r} = (med_gfp);

end
%%
