% function Nikon_Mig1Localization_26MAR2014()
% 
close all
clear all
close all hidden

clc

%Yoni
%folder_name_base = 'C:\Users\Yoni\Dropbox\gal_paper\Mig1localiztion\26MAR Mig1 localization DG\'

% folder_name_base='/Users/rae10/Dropbox/Galactose_Pathway/gal_paper/Mig1localiztion/26MAR Mig1 localization DG/'

load seq_cells_1_9
seg_data = seg_cells;
load seq_cells_10_11
seg_data{10} = seg_cells{10};
seg_data{11} = seg_cells{11};

load seq_cells_12_24
for i = 12:24
seg_data{i} = seg_cells{i};
end

for i = 1:24
    
%    folder_name = [folder_name_base,'\glucose_analog_C2_',num2str(i-1),'\']
    
     folder_name = [folder_name_base,'glucose_analog_C2_',num2str(i-1),'/']
    
    
    %files_gfp = dir([folder_name,'\*gfp*.tif']);
    
    files_gfp = dir([folder_name,'*gfp*.tif']);
    
    [val_bf,ind_bf] = sort([files_gfp.datenum]);
    files_gfp = files_gfp(ind_bf);
        
    
    I(:,:,i)  = mat2gray((imread([folder_name,files_gfp(1).name])));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     Construct a questdlg with three options
%     flag = 1;
%     temp = I(:,:,i);
%     k = 1;
%     while    flag==1
%         figure(1)
%         seg_cells{i}(:,:,k) = roipoly(temp);
%         
%         choice = questdlg('More cells?', ...
%             'Option Menu', ...
%             'Yes','No','Yes');
% 
%         switch choice
%             case 'Yes'
%                 flag = 1;
%                 temp = temp.*(~seg_cells{i}(:,:,k));
%             case 'No'
%                 flag = 0;
%         end
%         k = k+1;
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(i)
%     imshow(I(:,:,i));
%             choice = questdlg('More cells?', ...
%             'Continue?', ...
%             'Yes','No','Yes');
    
    II = mat2gray(I(:,:,i));
    bw = im2bw(II,graythresh(II));
    val_cell= II(find(bw));
    bw2 = im2bw(II,graythresh(val_cell));
%     figure;imshowpair(bw2,bw);
    I_proc(:,:,:,i) = imfuse(bw2,bw);
    val_cell = II(find(bw));
    cell(i) = length(val_cell);
    val_nuc = II(find(bw2));
    nuc(i) = length(val_nuc);
    ratio(i) = nuc(i)/cell(i);
    
    [t1,t2,num_cells] = size(seg_data{i});
    
    for j = 1:num_cells
        
        ratio_seg{i}(j) = length(find((bw2.*seg_data{i}(:,:,j))))/length(find((bw.*seg_data{i}(:,:,j))));
    end
    
end
%%

ind_m = [1 2 3 4 5 6;12 11 10 9 8 7;13 14 15 16 17 18;24 23 22 21 20 19]

% images

i = 1;I1=[];
for j = 1:6
    
    I1 = [I1 I(:,:,ind_m(i,j))];
end

i = 2;I2=[];
for j = 1:6    
    I2 = [I2 I(:,:,ind_m(i,j))];
end

i=3;I3=[];
for j = 1:6
    
    I3 = [I3 I(:,:,ind_m(i,j))];
end

i = 4;I4=[];
for j = 1:6
    
    I4 = [I4 I(:,:,ind_m(i,j))];
end

I_final = [I1;I2;I3;I4];
% proceced imges

i = 1;I1=[];
for j = 1:6
    
    I1 = [I1 I_proc(:,:,:,ind_m(i,j))];
end

i = 2;I2=[];
for j = 1:6    
    I2 = [I2 I_proc(:,:,:,ind_m(i,j))];
end

i=3;I3=[];
for j = 1:6
    
    I3 = [I3 I_proc(:,:,:,ind_m(i,j))];
end

i = 4;I4=[];
for j = 1:6
    
    I4 = [I4 I_proc(:,:,:,ind_m(i,j))];
end
I_final_proc = [I1;I2;I3;I4];


figure(1)
imshow(I_final,[]);
figure(2)
imshow(I_final_proc,[]);

%%

%figure_folder = 'C:\Users\ys151\Dropbox\gal_paper\figures\Matlab Source\';
figure_folder = '.'

cmap=cbrewer('seq', 'Blues', 25);

ratio_M = [ratio(1:6);fliplr(ratio(7:12));ratio(13:18);fliplr(ratio(19:24))];

ratio_M = 1- mat2gray(ratio_M);
figure(1);plot(ratio_M','o-','linewidth',2);xlim([0.5,6.5]);
set(gca,'xtick',([1:6]),'xticklabel',['No',  mat2cell(-7:2:1)]);
Set_fig_YS(figure(1),18,12,12);
print('-dtiff','-r1000',[figure_folder,'mig1localizationraw']);


ratio_M = padarray(ratio_M,[1 1]);
figure(3)
h = pcolor(ratio_M(2:end,2:end));colormap(cmap);colorbar
% ylabel('Glucose [log_{2}(%)] ')
% xlabel('Galactose [log_{2}(%)] ')
% title('Mig1p localization','fontsize',22)

set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 4:1:25]+0.5,'xticklabel',['No', mat2cell(-7:2:1)]);
set(gca,'ytick',[1 2 3 4:1:22]+0.5,'yticklabel',['No', mat2cell(-7:2:-3)]);

colorbar 
colormap(cmap)

Set_fig_YS(figure(3),18,18,22)
filename=['Mig1p_localization'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
%%


cmap=cbrewer('seq', 'Blues', 25);

ratio_count  = [1
0.923076923
0.705882353
0.785714286
0.684210526
0.891891892
0.852941176
0.8
0.666666667
0.714285714
0.64516129
0.709677419
0.095238095
0.166666667
0.076923077
0.3
0.090909091
0.068965517
0.1
0
0
0
0
0
]';


ratio_M_count = [ratio_count(1:6);fliplr(ratio_count(7:12));ratio_count(13:18);fliplr(ratio_count(19:24))];
figure(4);plot(ratio_M_count','o-','linewidth',2);xlim([0.5,6.5]);
set(gca,'xtick',([1:6]),'xticklabel',['No',  mat2cell(-7:2:1)]);
%Set_fig_YS(figure(4),18,12,12);
%print('-dtiff','-r1000',[figure_folder,'mig1localizationraw']);

% ratio_M = 1- mat2gray(ratio_M);
ratio_M_count = padarray(ratio_M_count,[1 1]);
figure(5)
h = pcolor(ratio_M_count(2:end,2:end));colormap(cmap);colorbar
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1 2 3 4:1:25]+0.5,'xticklabel',['No|',sprintf('%0.1e|', (2.^(-7:2:1)))]);
set(gca,'ytick',[1 2 3 4:1:22]+0.5,'yticklabel',['No|', sprintf('%0.1e|',2.^(-7:2:-3))]);

colorbar 
colormap(cmap)
Set_fig_RE(figure(5),18,18,22)

% filename=['Mig1p_localization.pdf']
% export_fig(filename, '-pdf','-transparent','-nocrop');

%%
average_localization=sum(ratio_M_count,2)./6;

average_localization=[average_localization average_localization];

average_localization=average_localization(2:end,:);
hfig=figure('Position',[  440   378    92   420]);
h=pcolor(average_localization);
colormap(cmap);
%colorbar
set(h,'edgecolor','none');
Set_fig_RE(hfig,18,18,22)

filename=['Mig1p_localization_average.pdf']
export_fig(filename, '-pdf','-transparent','-nocrop');


%%
th = 0.4;
figure(10)
ind = [1:6,12:-1:7,13:18,24:-1:19]
for i = 1:24
    
    figure(10)
    subplot(4,6,ind(i))
    hist(ratio_seg{i});
    title(['m  = ',num2str(mean(ratio_seg{i})),' f= ',num2str(ratio(ind(i)))]);
    xlim([0 0.7]);
    figure(11);hold on;
    [y,x] = hist(ratio_seg{i});
    bar(x,y);
    num_cells = length(ratio_seg{i});
    ratio_fration(i) =  sum(ratio_seg{i}<th)/num_cells;
    
    
end

%%
% [Final1,rect] = imcrop(I_final);
close all

rect = [1.9595    0.6205    0.1800    0.1620]*1e3;
final1 = imcrop(I_final,rect);
final1_proc = imcrop(I_final_proc,rect)
figure;imshow(final1,[]);
figure;imshow(final1_proc);


rect2 = [1.8215    2.5945    0.1500    0.2040]*1e3;
final2 = imcrop(I_final,rect2);
final2_proc = imcrop(I_final_proc,rect2)
figure;imshow(final2,[]);
figure;imshow(final2_proc);

% end

% 