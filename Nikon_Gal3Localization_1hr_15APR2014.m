close all
clear all
close all hidden
clc


folder_name_base = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\Gal3Localization4H_15APRR2014\'
I1 = [];

for i = 1:6
    
    folder_name = [folder_name_base,'\glucose_analog_C2_',num2str(i-1),'\']
    
    
    
    files_gfp = dir([folder_name,'\*gfp*.tif']);
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);
    
    
    
    I1 = [I1 mat2gray((imread([folder_name,files_gfp(1).name])))];
    
    
end
I2 = [];

for i = 7:12
    
    folder_name = [folder_name_base,'\glucose_analog_C2_',num2str(i-1),'\']
    
    
    
    files_gfp = dir([folder_name,'\*gfp*.tif']);
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);
    
    
    
    I2 = [I2 mat2gray((imread([folder_name,files_gfp(1).name])))];
    
    
end
%%
figure
I = [I1;fliplr(I2)];
imshow(I)
% rep = 3
% start(1,1,:) = ones(1,3)*300;
% start(2,1,:) = ones(1,3)*400;
% start(3,1,:) = ones(1,3)*550;

% ratio_M = [ratio(1:6);fliplr(ratio(7:12));ratio(13:18);fliplr(ratio(19:24))]
% ratio_M = padarray(ratio_M,[1 1]);
