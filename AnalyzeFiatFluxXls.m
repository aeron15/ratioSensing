%% Read the excel output of FiatFlux

clear all
close all
clc

% Yoni
%folder_name = 'C:\Users\ys151\Dropbox\C13 data\140507_Springer\'

% Renan
folder_name='/Users/ys151/Dropbox/C13 data/140507_Springer/'

filesname = dir([folder_name,'*.xls'])


for i = 1:length(filesname)
[num,txt,raw] = xlsread([folder_name,filesname(i).name]);
aa(i,:) = num(62:end,1);
temp = aa(i,find(aa(i,:)));
temp = temp(find(~isnan(temp)));
aaa{i} = temp;
thresh(i) = graythresh(aaa{i});
ind = find(aaa{i}>thresh(i))
mode_aa(i) = median(aaa{i}(ind))
figure(i)
hist(aaa{i});hold on

end
%%

ratio_c13_13{1} =mode_aa(1:3);
ratio_c13_13{2} =mode_aa(4:5);
ratio_c13_13{3} =mode_aa(6:8);
ratio_c13_13{4} =mode_aa(9:11);
ratio_c13_13{5} =mode_aa(12);
ratio_c13_13{6} =mode_aa(13:15);
ratio_c13_13{7} =mode_aa(16:18);
ratio_c13_13{8} =mode_aa(19:21);

%%
glc = [-4 -6.5];gal = [-5 -1 1];
[GAL,GLC] = meshgrid(gal,glc);
GAL = 2.^GAL;
GLC = 2.^GLC;

data = [mode_aa(1), mode_aa(3), mode_aa(5);
        mode_aa(2), mode_aa(4), nan]
gal_incrop = 1-data;

plot(gal,gal_incrop)
% 1 3 5
% 2 4 6

x = 55*2.^gal(1:2)';y=55*2.^glc';
[XOut, YOut, ZOut] = prepareSurfaceData(x,y,gal_incrop(1:2,1:2));

s = fit([XOut YOut],ZOut,'(v/y)*1/(1 + kgal/x*(1+y/kglc))') 


% 
% 
% glc_ratio = gal_incrop(2,:)./gal_incrop(1,:)
% 
% data = [mode_aa(1), mode_aa(3), mode_aa(5);
%         mode_aa(2), mode_aa(4), nan]
% gal_incrop = 1-data;
% 

% 
%  
s = fit(2.^gal',gal_incrop(1,:)','v/(1 + kgal/x)')
s2 = fit(2.^gal',gal_incrop(1,:)','v/(1 + kgal/x*(1+0.0625/kglc))')


syms kgal kglc
kgal = 1
kglc = 0.05
 v = 1./ ( 1 + (kgal./GAL).*(1+GLC./kglc)   )
 v(1,:)./v(2,:)
%  
s = solve( (v(2,1)/v(1,1))/ (v(2,2)/v(1,2)) == (gal_incrop(2,1)/gal_incrop(1,1))/(gal_incrop(2,2)/gal_incrop(1,2)),kglc,'Real',true)
%  


% % 
% figure(1)
% plot(2.^gal,gal_incrop(1,:),'.-');hold on;
% plot(2.^gal,s(2.^gal),'r');
% plot(2.^gal,s2(2.^gal),'g');
% 
% s3 = fit(2.^gal',gal_incrop(1,:)','v/(1 + kgal/x*(1+0.062/kglc))')
% 
% 
% glc_ratio = gal_incrop(2,:)./gal_incrop(1,:);
% 
