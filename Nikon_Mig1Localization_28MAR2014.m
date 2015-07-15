close all
clear all
close all hidden
clc


folder_name_base = 'C:\Users\Savir\Desktop\Mig localizationDG\Mig localizationDG'
pos = [ 1 3 5 2 4 6]
for i = 1:24
    
    folder_name = [folder_name_base,'\glucose_analog_C2_',num2str(i-1),'\']
    
    
    
    files_gfp = dir([folder_name,'\*gfp*.tif']);
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);
    
    
    
    I(:,:,i) =  ((imread([folder_name,files_gfp(1).name])));
    
    II = mat2gray(I(:,:,i));%figure;imshow(II,[]);
    % [y,c] = kmeans(II(:),3,'replicates',3,'start',start);
    bw = im2bw(II,graythresh(II));%figure;imshowpair(bw,II)
    val_cell= II(find(bw));
    bw2 = im2bw(II,graythresh(val_cell));%figure;imshowpair(bw2,bw)
    val_cell = II(find(bw));
    cell(i) = length(val_cell);
    val_nuc = II(find(bw2));
    nuc(i) = length(val_nuc);
    ratio(i) = nuc(i)/cell(i);
    
    
    %     figure(i)
    %     subplot(3,2,pos(i))
    %     imshow(I(:,:,i),[])
    
    %     x = I(:,:,i);x = double(x(:));
    %     figure(2)
    %     subplot(3,2,pos(i))
    %     hist(x,10)
    %     [y,c] = kmeans(x,3,'replicates',3);
    %     [temp,ind] = sort(c);
    %     figure(3)
    %     subplot(2,3,i)
    %     bar([sum(y==ind(1)),sum(y==ind(2)),sum(y==ind(3))])
    
end
%%
% rep = 3
% start(1,1,:) = ones(1,3)*300;
% start(2,1,:) = ones(1,3)*400;
% start(3,1,:) = ones(1,3)*550;

ratio_M = [ratio(1:6);fliplr(ratio(7:12));ratio(13:18);fliplr(ratio(19:24))]
ratio_M = padarray(ratio_M,[1 1]);


