% close all
% clear all
% clc
% %%
% %folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\17APR Gal1 dynamics\images\'
% folder_name='/Volumes/sysbio/Springer Lab/Lab Members/Yoni/Nikon/17APR Gal1 dynamics/images/'
%     
%pos_num = 20;
dt=15/60;
% ind_position = [1:pos_num];
% %%
% close all
% for i = 1:length(ind_position)
%     
%     
%     files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'/*rfp*.tif']);
%     [val_bf,ind_bf] = sort([files_rfp.datenum]);
%     files_rfp = files_rfp(ind_bf);
%     
%     
%     files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'/*yfp*.tif']);
%     [val_bf,ind_bf] = sort([files_yfp.datenum]);
%     files_yfp = files_yfp(ind_bf);
%     
%     j0=0;
%     for j = 1:50
%         [i,j]
%         
%         
%         I_yfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'/',files_yfp(j+j0).name]));
%         
%         I_rfp =  (imread([folder_name,'pos',num2str(ind_position(i)),'/',files_rfp(j+j0).name]));
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
% % 
% % 
% % save('C:\Users\ys151\Desktop\WorkingZone\Gal\Gal1DynamicsOnix_17_APR_2014\data.mat')
% 
%%
%load('C:\Users\ys151\Desktop\WorkingZone\Gal\Gal1DynamicsOnix_17_APR_2014\data.mat')
load('Gal1DynamicsOnix_17_APR_2014/data.mat')

%%
close all
figure(1)

T =[0:49]*dt
plot(T,area_rfp,'.-r');ylim([0 1e6])

%% Calcualte growth rate

x = [20:40]; Y = area_rfp(:,21:41);
ind_pos = [1:20];
for i = 1:length(ind_pos)
    s = fit(x',log(Y(i,:))','a*x+b','robust','on');
    m(i) = s.a;
    C(:,:,i) = confint(s);
    
     figure(5)
    subplot(4,5,i);
    plot(x,log(Y(i,:))','.');hold on;
    plot(x,s(x));
    
end
figure(2)
eu = squeeze(C(2,1,:))';
ed = squeeze(C(1,1,:))';
errorbar(ind_pos,m,m-ed,eu-m,'o');

    


%%
close all
ctrs = linspace(80,1000,200);

clear yfp_hist_conc med rfp_hist_conc med_r yfp_norm_hist_conc med_norm
for i = 1:pos_num
    for j=1:50
        
        yfp_hist_conc(i,j,:) = hist(yfp_hist{i,j},ctrs);
        yfp_hist_conc(i,j,:) = yfp_hist_conc(i,j,:) /sum(yfp_hist_conc(i,j,:) );
        med(i,j) = median(yfp_hist{i,j});
        
        rfp_hist_conc(i,j,:) = hist(rfp_hist{i,j},ctrs);
        rfp_hist_conc(i,j,:) = rfp_hist_conc(i,j,:) /sum(rfp_hist_conc(i,j,:) );
        med_r(i,j) = median(rfp_hist{i,j});
        
        yfp_norm_hist_conc(i,j,:) = hist(yfp_hist{i,j});
        yfp_norm_hist_conc(i,j,:) = yfp_norm_hist_conc(i,j,:) /sum(yfp_norm_hist_conc(i,j,:) );
        med_norm(i,j) = median(yfp_hist{i,j});
        
        
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

time_probe = [1 2 3 4 16 32]


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
    for j=1:50
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
time_probe = [0 1 2 4 8 16 32 40]+1;

TT = T(time_probe);

k =1;
for i = 1:4
    x = squeeze(yfp_hist_conc(i,time_probe,:));
    
    figure(3)
    for j = 1:length(TT)
        subplot(4,length(time_probe),k)
        %plot(ctrs,x(1,:),'r','linewidth',2);hold on;
        plot(ctrs,x(j,:),'k','linewidth',2);set(gca,'xscale','log');xlim([80,1000]);ylim([0 1.1]);title(num2str((TT(j))));
        axis square
        box off
        set(gcf,'color','None')
        set(gca,'YTick',[],'YTickLabel','')
    
            k = k+1;

    end
end
    %Set_fig_YS(figure(3),12,18,18)

%% Print images for micrographs

% close all

time_probe = [4 8 16 32 40]+1;

%time_probe=32;
TT = T(time_probe);
% 5 7 14 19
ind_position = [5 7 14 19];
%ind_position = [5]

I_yfp=[];
I_rfp=[];
I_bf=[];
II_rfp=[];
II_yfp = [];
II_bf=[];
for i = 1:length(ind_position)
    
    %folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\17APR Gal1 dynamics\images\';
    folder_name='/Volumes/sysbio/Springer Lab/Lab Members/Yoni/Nikon/17APR Gal1 dynamics/images/'

    files_bf = dir([folder_name,'pos',num2str(ind_position(i)),'/*bf*.tif']);
    [val_bf,ind_bf] = sort([files_bf.datenum]);
    files_bf = files_bf(ind_bf);
    
    
    files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'/*yfp*.tif']);
    [val_bf,ind_bf] = sort([files_yfp.datenum]);
    files_yfp = files_yfp(ind_bf);
    
    files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'/*rfp*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
    
    
    j0=0;
    I_yfp=[];
    I_rfp=[];
    I_bf=[];
    for j = 1:length(time_probe)
        [i,j]
        
        
        I_yfp = [I_yfp, (imread([folder_name,'pos',num2str(ind_position(i)),'/',files_yfp(time_probe(j)).name]))];
        
        I_rfp = [I_rfp, (imread([folder_name,'pos',num2str(ind_position(i)),'/',files_rfp(time_probe(j)).name]))];
        
        I_bf = [I_bf, mat2gray(imread([folder_name,'pos',num2str(ind_position(i)),'/',files_bf(time_probe(j)).name]))];

    end
        
    II_bf = [II_bf;I_bf];
    II_rfp = [II_rfp;I_rfp];
    II_yfp = [II_yfp;I_yfp]; 
%     I_yfp = mat2gray(I_yfp);

     
    
    %figure
%     imshowpair(I_yfp,I_rfp,'method','blend','scaling','none');
    %imshow(I_rfp,[0 1]);

end

%%

    II_yfp = mat2gray(II_yfp);
    II_rfp = mat2gray(II_rfp);

    III_bf = 0* repmat(II_bf,[1 1 3]);
    
    III_bf(:,:,2) = III_bf(:,:,2)+II_yfp;
    III_bf(:,:,1) = III_bf(:,:,1)+II_rfp;
    
    figure(4)
    imshow(III_bf);
    
%% 20140701 Generate figure 1B for Mike

close all;
%time_probe = [0 1 2 4 8 16 32 40]+1;
time_probe = [32]+1;

TT = T(time_probe);

k =1;
for i = 1:3
    x = squeeze(yfp_hist_conc(i,time_probe,:));
    
    x=yfp_hist_conc(i,time_probe,:);
    
    figure(3)
    for j = 1:length(TT)
        subplot(4,length(time_probe),k)
        %plot(ctrs,x(1,:),'r','linewidth',2);hold on;
        plot(ctrs,x(j,:),'k','linewidth',2);
        set(gca,'xscale','log');xlim([80,1000]);ylim([0 1.1]);
        %title(num2str((TT(j))));
        axis square
        box off
        set(gcf,'color','None')
        set(gca,'YTick',[],'YTickLabel','')
    
            k = k+1;

    end
end

filename='output/Figure_1B_microscopy_pannels';
export_fig(filename, '-pdf','-transparent','-nocrop');
