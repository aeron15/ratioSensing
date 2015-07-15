%% Analyze time dynamics of H3 and A10 response


close all
clear all
clc


%%
load 'C:\Users\ys151\Desktop\WorkingZone\Gal\FACS time points\screen';
well_num_str= {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12' };
well_letter = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};

%%
ctrs = logspace(-4,6,200);
names = fieldnames(screen);
A10_ind = [1:7];
H3_ind = [8:14];
data = struct2cell(screen);

for i = 1:length(A10_ind)
    wells = fieldnames(data{i});
    for j = 1:length(wells)
        x = cell2mat(wells(1));
        r = find(strcmp(x(1),well_letter));
        c = find(strcmp(x(2:3),well_num_str));
        figure(i)
        subplot(8,12,c+(r-1)*12);
        [y,x] = hist(
