close all
clear all
clc

folder_name = 'C:\Users\ys151\Desktop\16_1_2013\'
pos_num = 40;
dt=20/60;

%%
close all

for i = 1:pos_num
    tic
    files = dir([folder_name,'*r FITC_s',num2str(i),'_*.tif']);
    [val,ind] = sort([files.datenum]);
    files = files(ind);
    
    files_bf = dir([folder_name,'*e FITC_s',num2str(i),'_*.tif']);
    [val_bf,ind_bf] = sort([files_bf.datenum]);
    files_bf = files_bf(ind_bf);
    
    
    files_rfp = dir([folder_name,'*r TRITC_s',num2str(i),'_*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
    files_bfp = dir([folder_name,'*r DAPI_s',num2str(i),'_*.tif']);
    [val_bf,ind_bf] = sort([files_bfp.datenum]);
    files_bfp = files_bfp(ind_bf);
    j0=11;
    for j = 1:length(files)-j0
        [i,j]
        I_bf(:,:,j) = (imread([folder_name,files_bf(j+j0).name]));
        
        
        I_yfp(:,:,j) =  (imread([folder_name,files(j+j0).name]));
        
        I_rfp(:,:,j) = (imread([folder_name,files_rfp(j+j0).name]));
        I_bfp(:,:,j) = (imread([folder_name,files_bfp(j+j0).name]));
        
        
        
    end
    
    %   I_bf = mat2gray(I_bf);I_yfp = mat2gray(I_yfp);I_rfp = mat2gray(I_rfp);I_bfp = mat2gray(I_bfp);
   [k,l,t] = size(I_yfp)
   
   h = fspecial('average',3);
   I = imfilter(I_bfp(:,:,1),h);
   graythresh(mat2gray(I))

b=140;
    for j = 1:t
        [i,j,2]
        BW_rfp = im2bw( mat2gray(I_rfp(:,:,j)-b) ,graythresh(mat2gray(I_rfp(:,:,j))-b ));
        BW_yfp = im2bw( mat2gray(I_yfp(:,:,j) ),graythresh(mat2gray(I_yfp(:,:,j))));
        BW_bfp = im2bw(mat2gray(I_bfp(:,:,j)) ,graythresh(mat2gray(I_bfp(:,:,j))));%graythresh(mat2gray(I_bfp(:,:,j),[b,double(max(max(I_bfp(:,:,j))))])));
        r_bfp(i,j) = sum(BW_bfp.*BW_yfp)/sum(BW_bfp);
        r_rfp(i,j) = sum(BW_rfp.*BW_yfp)/sum(BW_rfp);
        area_bfp(j) = length(find(BW_bfp));
        area_rfp(j) = length(find((BW_rfp)));
%        [r,c] = find(BW_yfp);
%        int_yfp(j) = sum(sum(I_yfp(r,c,j)));
    end
    %     [i,j,2]
    %     x(:,:,1) = I_bf(:,:,j);x(:,:,2) = I_bf(:,:,j);x(:,:,3) = I_bf(:,:,j);
    %     y(:,:,2) = I_yfp(:,:,j);y(:,:,1) = 0;y(:,:,3)=0;
    %     z(:,:,1) = I_rfp(:,:,j);z(:,:,2) = 0;z(:,:,3)=0;
    %     w(:,:,3) = I_bfp(:,:,j);w(:,:,2) = 0;w(:,:,2)=0;
    %
    %
    % end
%      M(:,:,1,:) = mat2gray(I_rfp);
%         M(:,:,2,:) = mat2gray(I_yfp);
%         M(:,:,3,:) = mat2gray(I_bfp);
%     save(['16_1_2013_M_',num2str(i)],'M');
    
    toc
end
i
save 16_1_2013 area_bfp area_rfp r_bfp r_rfp
%%
close all
T = [0:t-1]*dt+j0*dt;
subplot(2,1,1)
plot(T,area_bfp,'b',T,area_rfp,'r');
subplot(2,1,2);
plot(r_bfp','b');hold on;plot(r_rfp','r');

%%
ctrs = [0:0.05:1];
[c1,x1] = hist(r_rfp(1:10,14),ctrs);
[c2,x2] = hist(r_rfp(11:20,14),ctrs);
[c3,x3] = hist(r_rfp(21:30,14),ctrs);
[c4,x4] = hist(r_rfp(31:40,14),ctrs);

figure(2);hold on;
bar(x1,c1,'r');
bar(x2,c2,'g');
bar(x3,c3,'b');
bar(x4,c4,'k');

 median(r_rfp(1:10,14))
 median(r_rfp(11:20,14))
median(r_rfp(21:30,14))
median(r_rfp(31:40,14))
