close all
clear all
clc
%%
% folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\20140409_Gal1_dynamics\09APR Gal1 induction\images\'
% pos_num = 20;
% dt=15/60;
% ind_position = [1:pos_num];
% %%
% close all
% for i = 1:length(ind_position)
%     
%     
%     files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*rfp*.tif']);
%     [val_bf,ind_bf] = sort([files_rfp.datenum]);
%     files_rfp = files_rfp(ind_bf);
%     
%     
%     files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
%     [val_bf,ind_bf] = sort([files_yfp.datenum]);
%     files_yfp = files_yfp(ind_bf);
%     
%     j0=0;
%     for j = 1:length(files_rfp)-j0
%         [i,j]
%         
%         
%         I_yfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(j+j0).name]));
%         
%         I_rfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(j+j0).name]));
%         
%         BW_rfp = im2bw( mat2gray(I_rfp) ,graythresh(mat2gray(I_rfp )));
%         BW_yfp = im2bw( mat2gray(I_yfp),graythresh(mat2gray(I_yfp)));
%         
%         area_yfp(i,j) = length(find(BW_yfp));
%         area_rfp(i,j) = length(find((BW_rfp)));
%         %yfp
%         xx = double(I_yfp);
%         yy = xx.*BW_rfp;
%         ind = find(yy>0);        
%         yfp_hist{i,j} = yy(ind);
%         
%         %rfp
%         xx = double(I_rfp);
%         yy = xx.*BW_rfp;
%         ind = find(yy>0);
%         rfp_hist{i,j} = yy(ind);
%         
%     end
%     
% end
% 
% 
% 
% k = length(area_rfp);
% folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\20140409_Gal1_dynamics\09APR Gal1 induction2\images\'
% 
% pos_num = 20;
% ind_position = [1:pos_num];
% 
% for i = 1:length(ind_position)
%     
%     
%     
%     files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*rfp*.tif']);
%     [val_bf,ind_bf] = sort([files_rfp.datenum]);
%     files_rfp = files_rfp(ind_bf);
%     
%     
%     files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
%     [val_bf,ind_bf] = sort([files_yfp.datenum]);
%     files_yfp = files_yfp(ind_bf);
%     
%     j0=0;
%     for j = 1:length(files_rfp)-j0
%         [i,j]
%         
%         
%         I_yfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(j+j0).name]));
%         
%         I_rfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(j+j0).name]));
%         
%         BW_rfp = im2bw( mat2gray(I_rfp) ,graythresh(mat2gray(I_rfp )));
%         BW_yfp = im2bw( mat2gray(I_yfp),graythresh(mat2gray(I_yfp)));
%         
%         area_yfp(i,k+j) = length(find(BW_yfp));
%         area_rfp(i,k+j) = length(find((BW_rfp)));
%         
%         %yfp
%         xx = double(I_yfp);
%         yy = xx.*BW_rfp;
%         ind = find(yy>0);
%         yfp_hist{i,k+j} = yy(ind);
%         %rfp
%         xx = double(I_rfp);
%         yy = xx.*BW_rfp;
%         ind = find(yy>0);
%         rfp_hist{i,k+j} = yy(ind);
%         
%     end
%     
% end

% save('C:\Users\ys151\Desktop\WorkingZone\Gal\Gal1DynamicsOnix_09_APR_2014\data.mat')

%%

%load('C:\Users\Tal Aizenfeld\Desktop\Gal1DynamicsOnix_09_APR_2014\data.mat')
close all
figure(1)

cutoff = 48;
dT = 31/60;
dt=0.25;

T1 =[0:23]*dt
T2 = [0:23]*dt+dT+T1(end);
T = [T1 T2];
plot(T,area_rfp(:,1:48),'.-r');ylim([0 1e6])
%%
close all
ctrs = linspace(80,1000,200);
clear yfp_hist_conc med rfp_hist_conc med_r yfp_norm_hist_conc med_norm
for i = 1:pos_num
    for j=1:48
        
        yfp_hist_conc(i,j,:) = hist(yfp_hist{i,j},ctrs);
        yfp_hist_conc(i,j,:) = yfp_hist_conc(i,j,:) /sum(yfp_hist_conc(i,j,:) );
        med(i,j) = median(yfp_hist{i,j});
        
        rfp_hist_conc(i,j,:) = hist(rfp_hist{i,j},ctrs);
        rfp_hist_conc(i,j,:) = rfp_hist_conc(i,j,:) /sum(rfp_hist_conc(i,j,:) );
        med_r(i,j) = median(rfp_hist{i,j});
        
        yfp_norm_hist_conc(i,j,:) = hist(yfp_hist{i,j}./rfp_hist{i,j},ctrs);
        yfp_norm_hist_conc(i,j,:) = yfp_norm_hist_conc(i,j,:) /sum(yfp_norm_hist_conc(i,j,:) );
        med_norm(i,j) = median(yfp_hist{i,j}./rfp_hist{i,j});
        
        
    end
    figure(2)
    subplot(4,5,i)
    x = squeeze(yfp_hist_conc(i,:,:));
    pcolor(T,ctrs,x');
    figure(3)
    subplot(4,5,i)
    x = squeeze(rfp_hist_conc(i,:,:));
    pcolor(T,ctrs,x');
    
    figure(4)
    subplot(4,5,i)
    plot(T,med(i,:));title(num2str(i));
    
    figure(5)
    subplot(4,5,i)
    plot(T,med_r(i,:));title(num2str(i));
    
    figure(6)
    subplot(4,5,i)
    plot(T,med_norm(i,:));title(num2str(i));
end

time_probe = [1 2 4 8 16 32 18]


for i = 1:pos_num
    x = squeeze(yfp_hist_conc(i,time_probe,:));
    
    figure(7)
    subplot(4,5,i)
    plot(ctrs,x')
    
    x = squeeze(rfp_hist_conc(i,time_probe,:));   
    figure(8)
    subplot(4,5,i)
    plot(ctrs,x')
    
end

%% close all
close all
ctrs = linspace(80,1000,100);
clear yfp_hist_conc med rfp_hist_conc med_r
for i = 1:4
    for j=1:48
        y = yfp_hist{5*(i-1)+1:5*(i-1)+5,j};
        yfp_hist_conc(i,j,:) = hist(y,ctrs);
        yfp_hist_conc(i,j,:) = yfp_hist_conc(i,j,:) /max(yfp_hist_conc(i,j,:) );
        med(i,j) = median(y);
        
        y = rfp_hist{5*(i-1)+1:5*(i-1)+5,j};
        rfp_hist_conc(i,j,:) = hist(y,ctrs);
        rfp_hist_conc(i,j,:) = rfp_hist_conc(i,j,:) /max(rfp_hist_conc(i,j,:) );
        med_r(i,j) = median(y);
        
    end
    figure(2)
    subplot(4,1,i)
    x = squeeze(yfp_hist_conc(i,:,:));
    h=pcolor(T,ctrs,x');set(h,'edgecolor','none');set(gca,'yscale','log');ylim([80,1000])
    figure(3)
    plot(T,med(i,:));hold on;set(gca,'yscale','linear');
    
end

%%

time_probe1 = [0 1 2 4 8 16]+1;
T1(time_probe1)
time_probe2 = [8 24];
T2(time_probe2)

time_probe = [time_probe1 time_probe2+24];

T(time_probe);
TT = T(time_probe);TT(end-1)=8;TT(end)=12;

for i = 1:4
    x = squeeze(yfp_hist_conc(i,time_probe,:));
    
    figure(i)
    for j = 1:length(time_probe)
        subplot(1,length(time_probe),j)
        plot(ctrs,x(1,:),'r','linewidth',2);hold on;plot(ctrs,x(j,:),'k','linewidth',2);set(gca,'xscale','log');xlim([80,1000]);ylim([0 1.1]);title(num2str((TT(j))));
        axis square
    end
    Set_fig_RE(figure(i),12,18,18)
end

%% Print images for micrographs

close all

time_probe1 = 17;%[1 2 4 8 16]+1;
T1(time_probe1)
time_probe2 = [8 24];
T2(time_probe2)

ind_position = [2 6 12 17];
I_yfp=[];
I_rfp=[];
II_rfp=[];
II_yfp = [];
for i = 1:length(ind_position)
    
    folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\20140409_Gal1_dynamics\09APR Gal1 induction\images\';
    
    files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*bf*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
    
    files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
    [val_bf,ind_bf] = sort([files_yfp.datenum]);
    files_yfp = files_yfp(ind_bf);
    
    j0=0;
    I_yfp=[];
    I_rfp=[];
    for j = 1:length(time_probe1)
        [i,j]
        
        
        I_yfp = [I_yfp, (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(time_probe1(j)).name]))];
        
        I_rfp = [I_rfp, mat2gray(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(time_probe1(j)).name]))];
        
    end
    
    folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\20140409_Gal1_dynamics\09APR Gal1 induction2\images\';
    
    files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*bf*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
    
    files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
    [val_bf,ind_bf] = sort([files_yfp.datenum]);
    files_yfp = files_yfp(ind_bf);
    
    for j = 1:length(time_probe2)
        [i,j]
        
        
        I_yfp = [I_yfp, (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(time_probe2(j)).name]))];
        
        I_rfp = [I_rfp, mat2gray(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(time_probe2(j)).name]))];
        
    end
    
    
    
    
    II_rfp = [II_rfp;I_rfp];
    II_yfp = [II_yfp;II_yfp]; 
    I_yfp = mat2gray(I_yfp);

    
    
    
%     I_rfp(:,:,2) = I_rfp(:,:,2)+I_yfp;
%     I_rfp = mat2gray(I_rfp);
    
    
    figure
%     imshowpair(I_yfp,I_rfp,'method','blend','scaling','none');
    imshow(I_rfp,[0 1]);

end






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
