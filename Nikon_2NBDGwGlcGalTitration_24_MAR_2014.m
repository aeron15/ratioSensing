close all
clear all
close all hidden
clc


folder_name_base = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\2NBDG_titration_24_3_2014'
for i = 1:6
 
    folder_name = [folder_name_base,'\glucose_analog_',num2str(i),'\']
    
    %%
    
    files_gfp = dir([folder_name,'\*gfp*.tif']);
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);
    
    files_rfp = dir([folder_name,'\*rfp*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_dapi= files_rfp(ind_bf);
    
    
    gal10x(:,:,2) =  ((imread([folder_name,files_gfp(1).name])));
    gal10x(:,:,1) =  ((imread([folder_name,files_rfp(1).name])));
    gal10x(:,:,3) = 0;
    %%
    
    gal10x_N(:,:,1) = mat2gray(gal10x(:,:,1));
    gal10x_N(:,:,2) = mat2gray(gal10x(:,:,2));
    gal10x_N(:,:,3) = 0;
    
    %%
    figure(i)
    subplot(2,3,1)
    imshow(gal10x_N)
    subplot(2,3,2)
    imshow(gal10x_N(:,:,1))
    subplot(2,3,3)
    imshow(gal10x_N(:,:,2))
    %%
    
    bw = im2bw( gal10x_N(:,:,1),  graythresh( gal10x_N(:,:,1) ));
    subplot(2,3,4)
    imshow(bw)
    subplot(2,3,5)
    imshow(bw.*gal10x_N(:,:,2));
    subplot(2,3,6)
    imshow(~bw.*gal10x_N(:,:,2));
    %%subplot(3,1,)
    
    I = bw.*double(gal10x(:,:,2));
    signal{i} = ( I(find(I>0)));
    I = ~bw.*double(gal10x(:,:,2));
    back{i} = ( I(find(I>0)));
    
end
%%
c = [	0 14.25121951	27.82380952	53.11818182	116.86	194.7666667];

for i =1:6
    s_mean(i) = mean(signal{i});
    s_std(i) = std(signal{i});
    b_mean(i) = mean(back{i});
    b_std(i) = std(back{i});
end

figure(1)
errorbar(c,s_mean,s_std,'g');hold on;
errorbar(c,b_mean,b_std,'k')
figure(2)
plot(c,s_mean-b_mean,'o-k')

%%
