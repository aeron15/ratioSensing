% close all
% clear all
% clc
%%
% folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\17APR Gal1 dynamics\images\'
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
%     for j = 1:50
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
% save('C:\Users\ys151\Desktop\WorkingZone\Gal\Gal1DynamicsOnix_17_APR_2014\data.mat')

%%
%load('C:\Users\ys151\Desktop\WorkingZone\Gal\Gal1DynamicsOnix_17_APR_2014\data.mat')
close all
%%
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

time_probe = [0 4 8 16 32 48]+1;

TT = T(time_probe);
color_vec = [1 0 0;0 0 0;0 0 1;0 1/2 0;0 0 1/2;1 0.5 0];
%%
close all
k =1;
for i = 1:4
    
    figure(i)
    for j = 1:length(TT)
        x = squeeze(yfp_hist_conc(i,time_probe(j),:));
        x_0 = squeeze(yfp_hist_conc(i,1,:));
%         subplot(1,length(time_probe),k)
%         plot(ctrs,x_0,'r','linewidth',2);
        hold on;plot(ctrs,x,'color',color_vec(j,:),'linewidth',2);set(gca,'xscale','log');xlim([80,1000]);ylim([0 1.1]);title(num2str((TT(j))));hold on;
        axis square
            k = k+1;
    Set_fig_RE(figure(i),12,18,18)

    end
end

%% Print images for micrographs

% close all

time_probe = [0 4 8 16 32 40 48]+1;

time_probe=[0 32 48]+1;
TT = T(time_probe);
% 5 7 14 19
ind_position = [5  14 ];
% ind_position = 19
I_yfp=[];
I_rfp=[];
I_bf=[];
II_rfp=[];
II_yfp = [];
II_bf=[];
rect = [700.5100  301.5100  392.9800  427.9800];
for i = 1:length(ind_position)
    x
    folder_name = '\\research.files.med.harvard.edu\sysbio\Springer Lab\Lab Members\Yoni\Nikon\17APR Gal1 dynamics\images\';
    
    files_bf = dir([folder_name,'pos',num2str(ind_position(i)),'\*bf*.tif']);
    [val_bf,ind_bf] = sort([files_bf.datenum]);
    files_bf = files_bf(ind_bf);
    
    
    files_yfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*yfp*.tif']);
    [val_bf,ind_bf] = sort([files_yfp.datenum]);
    files_yfp = files_yfp(ind_bf);
    
    files_rfp = dir([folder_name,'pos',num2str(ind_position(i)),'\*rfp*.tif']);
    [val_bf,ind_bf] = sort([files_rfp.datenum]);
    files_rfp = files_rfp(ind_bf);
    
    
    
    j0=0;
    I_yfp=[];
    I_rfp=[];
    I_bf=[];
    for j = 1:length(time_probe)
        [i,j]
        if i == 5
            I = (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(time_probe(j)).name]));
        end
        I_yfp = [I_yfp, (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(time_probe(j)).name]))];
        
        I_rfp = [I_rfp, (imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(time_probe(j)).name]))];
        
        I_bf = [I_bf, (mat2gray(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_bf(time_probe(j)).name])))];

        
                
%         I_yfp = [I_yfp, imcrop(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_yfp(time_probe(j)).name]),rect)];
%         
%         I_rfp = [I_rfp, imcrop(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_rfp(time_probe(j)).name]),rect)];
%         
%         I_bf = [I_bf, mat2gray(imcrop(imread([folder_name,'pos',num2str(ind_position(i)),'\',files_bf(time_probe(j)).name]),rect))];
% 
%         
    end
    
    II_bf = [II_bf;I_bf];
    II_rfp = [II_rfp;I_rfp];
    II_yfp = [II_yfp;I_yfp]; 
%     I_yfp = mat2gray(I_yfp);

%      
%     
%     figure
% %     imshowpair(I_yfp,I_rfp,'method','blend','scaling','none');
%     imshow(I_rfp,[0 1]);

end

    II_yfp = mat2gray(II_yfp);
    II_rfp = mat2gray(II_rfp);

    III_bf = 0* repmat(II_bf,[1 1 3]);
    
    III_bf(:,:,2) = III_bf(:,:,2)+II_yfp;
    III_bf(:,:,1) = III_bf(:,:,1)+II_rfp;

    figure(4)
    imshow(III_bf);
%% Precent cell induction as a funciton of time - based on area

for i = 1:4
    for j = 1:length(T)
    x = squeeze(yfp_hist_conc(i,j,:));
    x_0 = squeeze(yfp_hist_conc(i,1,:));
    [perc_const,perc_adaptive,perc_area_temp(i,j)] = AnalyzeHist_YS(x,ctrs,2.5,1,x_0,ctrs);


    end
end

    
figure(1)
plot(T,perc_area_temp','.-','linewidth',2)
xlim([0 8]);box off

Set_fig_RE(figure(1),18,12,12);
% figure_folder = 'C:\Users\ys151\Dropbox\gal_paper\figures\Matlab Source\';
% print('-dtiff','-r1000',[figure_folder,'gal1dynamicsonix']);

%% Pull histograms for figure 1
close all

t_well = 8; % time of mesurment
ind = find(T==t_well);

x = squeeze(yfp_hist_conc(1,ind,:));

figure(1)
plot(ctrs,x/sum(x),'k','linewidth',6);box off;xlim([0 1000]);set(gca,'xscale','log')

Set_fig_YS(figure(1),25,18,18);

x = squeeze(yfp_hist_conc(2,ind,:));
figure(2)
plot(ctrs,x/sum(x),'k','linewidth',6);box off;xlim([0 1000]);set(gca,'xscale','log')
Set_fig_YS(figure(2),25,18,18)

    
x = squeeze(yfp_hist_conc(3,ind,:));
figure(3)
plot(ctrs,x/sum(x),'k','linewidth',6);box off;xlim([0 1000]);set(gca,'xscale','log')
Set_fig_YS(figure(3),25,18,18)

