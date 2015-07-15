close all
clear all
close all hidden
clc

folder_name = 'N:\Springer Lab\Lab Members\Yoni\Nikon\20MAR Mig1 Localization3\images\'
pos_num = 17;
dt=2/60;
ind_position = [4]
%%

for i = 1:length(ind_position)
    tic
    
    files_gfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*gfp*.tif']);
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);

    
        
    files_dapi= dir([folder_name,'pos',num2str(ind_position(i)),'\*dapi*.tif']);
    [val_bf,ind_bf] = sort([files_dapi.datenum]);
    files_dapi= files_dapi(ind_bf);
    
    
    j0=0;
    for j = 1:length(files_gfp)-j0
        [i,j]
        
        
        I_gfp(:,:,j) =  mat2gray((imread([folder_name,'pos',num2str(ind_position(i)),'\',files_gfp(j+j0).name])));
        
                
        I_dapi(:,:,j) =  mat2gray((imread([folder_name,'pos',num2str(ind_position(i)),'\',files_dapi(j+j0).name])));
                
                
        
        
    end
end

implay(I_gfp)
implay(I_dapi)
