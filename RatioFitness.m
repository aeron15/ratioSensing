%% ratio_fitness

close all
clear all
clc

load Growth_Rate % A10 A10 H3 H3
load final_yield_plate_1_A10

glc = 2.^[-5:1];
glc = [0 glc]
gal = 2.^[-8:1];
gal = [0 gal]

for i =1:4
    fs1_e(i,:) = growth_rate(1:8,1,i);
    fs2_e(i,:) = growth_rate(1,1:11,i);
    
end

figure(1)
subplot(2,1,1);hold on;
plot(log2(glc),fs1_e(1:2,:),'.-k');title('A10')
plot(log2(gal),fs2_e(1:2,:),'.-r');

subplot(2,1,2);hold on;
plot(log2(glc),fs1_e(3:4,:),'.-k');title('H3')
plot(log2(gal),fs2_e(3:4,:),'.-r');

%% 
close all
glc_vec = logspace(-3,0.6,100);
gal_vec = logspace(-3,0.6,100);

f1 = 0.07*ones(size(glc_vec));
f2 = 0.07*ones(size(gal_vec));

ind = find(gal_vec>(2^-6));
f2(ind) = 0.07+0.28*(1./(1 + (2^-6.5)./(gal_vec(ind)-(2^-6)).^1));

ind = find(glc_vec>(2^-7));
f1(ind) = 0.07+0.45*(1./(1 + (2^-6.5)./(glc_vec(ind)-(2^-7)).^1));

min_f1 = min(f1);
max_f1 = max(f1);
f1 = (f1-min_f1)/(max_f1-min_f1);
f2 = (f2-min_f1)/(max_f1-min_f1);

figure(1)
% plot(log2(glc),fs1_e(3,:),'.k');title('H3');hold on
plot(log2(glc_vec),f1,'k','linewidth',2);hold on

% plot(log2(gal),fs2_e(3,:),'.r');
plot(log2(gal_vec),f2,'color',[1 165/255 0],'linewidth',2);
xlim([-8,0])
box off
Set_fig_YS(figure(1),18,18,18)


delay = 0.37;
f11 = [0 f1(2:end)-delay];
% f11 = f1-0.17;
f11(f11<=0)=0;

figure(2)
plot(log2(glc_vec),f1,'k','linewidth',2);hold on;
plot(log2(gal_vec),f2,'color',[1 165/255 0],'linewidth',2);
plot(log2(glc_vec),f11,'--k','linewidth',2);hold on;
box off

% plot(log2(glc_vec(1:end)),f11,'--k')
xlim([-8,0])

Set_fig_YS(figure(2),18,18,18)

%%

cmap =  buildcmap_YS([0 0 0.5;1 1 1;0.5 0 0]);
colormap(cmap);
[s1,s2] = meshgrid(glc_vec,gal_vec);

for i = 1:length(glc_vec)
    for j = 1:length(gal_vec)
        
        if f1(i)>=f2(j)
            F(i,j) = f1(i);
            ind(i,j) = 0;
        else
            F(i,j) = f2(j);
            ind(i,j) = 1;
        end
    end
end

figure(3)
subplot(1,2,2)
h = surf(log2(gal_vec),log2(glc_vec),F);hold on;
set(h,'edgecolor','none','facecolor','interp');
view(0,90);axis square
xlim([-8 1]);ylim([-6 0]);
cmap =  buildcmap_YS([0 0 0.5;1 1 1;0.5 0 0]);

subplot(1,2,1)
h = surf(log2(gal_vec),log2(glc_vec),ind);hold on;
set(h,'edgecolor','none','facecolor','interp');
view(0,90);axis square
xlim([-8 1]);ylim([-6 0]);
colormap(cmap)
%


cost = 1
for i = 1:length(glc_vec)
    for j = 1:length(gal_vec)
        
        [val,ind_future] = min(abs(glc_vec-(glc_vec(i)/delay)));
        
        
        if f2(j) > f11(i)
            F(i,j) = f2(j);
            ind(i,j) = 1;
        else
            F(i,j) = f1(i);
            ind(i,j) = 0;
        end
%         if i <=delay
%             
%             if f2(j)>0.07
%                 F(i,j) = f2(j);
%                 ind(i,j) = 1;
%             else
%                 F(i,j) = 0.07;
%                 ind(i,j) = 0;
%             end
%         end  
%         if i>delay
%             
%             if f2(j)>f1(i-delay)
%                 F(i,j) = f2(j);
%                 ind(i,j) = 1;
%             else
%                 F(i,j) = f1(i);
%                 ind(i,j) = 0;
%             end
            
            
%         end
    end
end

figure(4)
subplot(1,2,2)
h = surf(log2(gal_vec),log2(glc_vec),F);hold on;
set(h,'edgecolor','none','facecolor','interp');
xlim([-8 1]);ylim([-6 0]);

view(0,90);axis square
subplot(1,2,1)
h = surf(log2(gal_vec),log2(glc_vec),ind);hold on;
set(h,'edgecolor','none','facecolor','interp');
view(0,90);axis square
xlim([-8 1]);ylim([-6 0]);
colormap(cmap);


Set_fig_YS(3,18,18,18);
Set_fig_YS(4,18,18,18);