close all
clear all
close all hidden
clc

for i = 1:3
    if i==1
        folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\2NBDG21MAR2014\GLC_10X\'
    end
    if i==2
        folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\2NBDG21MAR2014\No_GLC\'
    end
    if i==3
        folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\2NBDG21MAR2014\GAl_10X\'
    end
    
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
    figure(1)
    subplot(1,3,1)
    imshow(gal10x_N)
    subplot(1,3,2)
    imshow(gal10x_N(:,:,1))
    subplot(1,3,3)
    imshow(gal10x_N(:,:,2))
    %%
    
    bw = im2bw( gal10x_N(:,:,1),  graythresh( gal10x_N(:,:,1) ));
    figure(2)
    subplot(1,3,1)
    imshow(bw)
    subplot(1,3,2)
    imshow(bw.*gal10x_N(:,:,2));
    subplot(1,3,3)
    imshow(~bw.*gal10x_N(:,:,2));
    %%subplot(3,1,)
    
    figure(3)
    subplot(3,1,i)
    I = bw.*double(gal10x(:,:,2));
    [y,x] = hist( I(find(I>0)));hold on;
    bar(x,y/sum(y),'facecolor','g')
    I = ~bw.*double(gal10x(:,:,2));
    [y,x] = hist( I(find(I>0)));hold on;
    bar(x,y/sum(y),'facecolor','k');title('No GLC/GAL');xlim([120 260]);
    
end