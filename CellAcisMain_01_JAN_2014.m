close all
clear all
clc

folder_name = '\\N:\Springer Lab\Lab Members\Yoni\Nikon\20MAR Mig1 Localization3\images\pos4\'
pos_num = 17;
dt=20/60;
ind_position = [4]
%%
close all
ctrs = linspace(70,350,50);
for i = 1:length(ind_position)
    tic
    
    files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*rfp*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
   
    files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
    [val_bf,ind_bf] = sort([files_yfp.datenum]);
    files_yfp = files_yfp(ind_bf);
    
    j0=0;
    for j = 1:length(files_rfp)-j0
        [i,j]
        
        
        I_yfp(:,:,j) =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(j+j0).name]));
        
        I_rfp(:,:,j) =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(j+j0).name]));
        
        
        
    end
    
   [k,l,t] = size(I_yfp)
   


b=0;
    for j = 1:t
        BW_rfp = im2bw( mat2gray(I_rfp(:,:,j)-b) ,graythresh(mat2gray(I_rfp(:,:,j))-b ));
        BW_yfp = im2bw( mat2gray(I_yfp(:,:,j) ),graythresh(mat2gray(I_yfp(:,:,j))));

        area_yfp(i,j) = length(find(BW_yfp));
        area_rfp(i,j) = length(find((BW_rfp)));
        xx = double(I_yfp(:,:,j));
        yy = xx.*BW_rfp;
        ind = find(yy>0);
        yfp_hist(i,j,:) = hist(yy(ind),ctrs);
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
%%

ind_position = [2 8 10 11 16 17]

close all
figure(1)
T = [0:t-1]*dt+j0*dt;
plot(T(1:30),area_yfp(:,1:30),'k',T(1:30),area_rfp(:,1:30),'r');ylim([0 1e5])

for i = 1:length(ind_position)
    figure(i)
    x = squeeze(yfp_hist(i,:,:));
    pcolor(T(1:30),ctrs,x(1:30,:)')
end
    


%%
% ctrs = [0:0.05:1];
% [c1,x1] = hist(r_rfp(1:10,14),ctrs);
% [c2,x2] = hist(r_rfp(11:20,14),ctrs);
% [c3,x3] = hist(r_rfp(21:30,14),ctrs);
% [c4,x4] = hist(r_rfp(31:40,14),ctrs);
% 
% figure(2);hold on;
% bar(x1,c1,'r');
% bar(x2,c2,'g');
% bar(x3,c3,'b');
% bar(x4,c4,'k');
% 
%  median(r_rfp(1:10,14))
%  median(r_rfp(11:20,14))
% median(r_rfp(21:30,14))
% median(r_rfp(31:40,14))
