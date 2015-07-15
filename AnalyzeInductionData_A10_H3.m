function AnalyzeInductionData_A10_H3()

%% Analyze Deletion

load('../data/20140702_decision_front/output/data')
load('../data/20140702_decision_front/output/fit_17May')

fit_cuttoff = [2^-1 2^-7];
mid_value = 2^-4;
cutoff=0.1;

color_vec=cbrewer('qual', 'Dark2', 3);

color_vec=([0 0 1; 1 0 0])
%% Combine s288c H3 and A10
% close all
% 
% i=2
% WT{i}(find(WT{i}==-inf))=nan;
% WT{i}(find(WT{i}==0))=nan;
% 
% [x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT{i},max(max(WT{i})),min(min(WT{i})),cutoff,gal_s288b_4283other,glu_s288b_4283other,fit_cuttoff,mid_value);
% figure(2);
% x1=log2(x);
% y1=log2(y);
% idx_eliminate=y1>-0.1;
% y1(idx_eliminate)=[];
% x1(idx_eliminate)=[];
% plot(x1,y1,'.k','markersize',14);hold on;
% 
% 
% x2=x1;
% y2=s(x1);
% 
% idx_eliminate=y2<-7.1|y2>-1;
% y2(idx_eliminate)=[];
% x2(idx_eliminate)=[];
% 
% 
% plot(x2,y2,'k','linewidth',3);hold on;

%%
[x,y,s,a(2),b(2),a_d(2),a_u(2),b_d(2),b_u(2)] = SmoothHeatMap(A10,max(max(A10)),min(min(A10)),cutoff,gal,glu,fit_cuttoff,mid_value);

figure(2);
plot(log2(x),log2(y),'.','color',color_vec(1,:),'markersize',14);hold on;

% x1=log2(x);
% y1=s(log2(x));

x1=-8:0.1:0;
y1=s(x1);

idx_eliminate=y1<-7.1|y1>-1.1;
y1(idx_eliminate)=[];
x1(idx_eliminate)=[];

plot(x1,y1,'color',color_vec(1,:),'linewidth',3)

%%
[x,y,s,a(3),b(3),a_d(3),a_u(3),b_d(3),b_u(3)] = SmoothHeatMap(H3,max(max(H3)),min(min(H3)),cutoff,gal,glu,fit_cuttoff,mid_value);
figure(2);
plot(log2(x),log2(y),'.','color',color_vec(2,:),'markersize',14);
hold on;
x1=log2(x);
y1=s(log2(x));

idx_eliminate=y1<-7.1|y1>-1;
y1(idx_eliminate)=[];
x1(idx_eliminate)=[];


plot(x1,y1,'color',color_vec(2,:),'linewidth',3)


xlim([-8 0]);ylim([-8 0]);
axis square;
% legend({'s288c' 'A10' 'H3'})
ylabel('Glucose [log_{2}(%)] ')
xlabel('Galactose [log_{2}(%)] ')

Set_fig_RE(figure(2),18,18,22);

filename=['Decision_fronts_H3_A10'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%%

% plot(log2(x_May17),log2(y_May17),'o','markerfacecolor','none','color',[0 0 0],'markersize',2);hold on
% plot(log2(x_May17),s_May17(log2(x_May17)),'color',[0 0 0],'linewidth',2); hold on;
% xlim([-9 2]);ylim([-9 2]);

end
