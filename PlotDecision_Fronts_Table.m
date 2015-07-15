function PlotDecision_Fronts_Table()
% Plot previous concentrations of glucose and galactose assayed in the past

load('../data/20140702_decision_front/output/data')
load('../data/20140702_decision_front/output/fit_17May')

color_table=[
    
0.8203    0.4102    0.1172;
0    0.5000         0;
0.5000    0.5000    0.5000;
0    0.5000    0.5000;
1 0 0;
0 0 1;
];

Marker_Size=10;

fit_cuttoff = [2^-1 2^-7];
mid_value = 2^-4;
cutoff=0.1;

color_vec=cbrewer('qual', 'Dark2', 3);

%% Combine s288c H3 and A10
close all

i=2
WT{i}(find(WT{i}==-inf))=nan;
WT{i}(find(WT{i}==0))=nan;

[x,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(WT{i},max(max(WT{i})),min(min(WT{i})),cutoff,gal_s288b_4283other,glu_s288b_4283other,fit_cuttoff,mid_value);
%hfig2=figure(2);
hfig2=figure('Position',[321   250   600   419])
hold all;
x1=log2(x);
y1=log2(y);
idx_eliminate=y1>-0.1;
y1(idx_eliminate)=[];
x1(idx_eliminate)=[];
%plot(x1,y1,'.k','markersize',14);hold on;

x2=x1;
y2=s(x1);

idx_eliminate=y2<-7.1|y2>-1;
y2(idx_eliminate)=[];
x2(idx_eliminate)=[];


plot(x2,y2,'k','linewidth',1.5);hold on;


%% Concentrations Bennet experiments

Glc_bennet=[0.25];
Gal_bennet=[0.25];


plot(log2(Gal_bennet),log2(Glc_bennet),'o','MarkerSize',Marker_Size,'color',color_table(5,:),'MarkerFaceColor',color_table(5,:));


%% Concentrations Bennet experiments with 2% raffinose and 0.2% galactose in the supplementary material
%figure;hold all;
%X_ticks=[10^-4 2*10^-4 5*10^-4 10^-3 2*10^-3 5*10^-3 10^-2 2*10^-2 5*10^-2 10^-1  2*10^-1 4*10^-1 10^0 2*10^0];

X_ticks=[2.2,
1.1,
0.55,
0.275,
0.1375,
0.06875,
0.034375,
0.0171875,
0.00859375,
0.004296875,
0.0021484375,
0.00107421875,
0.000537109375,
0.0002685546875,
0.0001342773438];

Y_ticks=repmat(0.00001,1,length(X_ticks));

plot(log2(X_ticks),log2(Y_ticks),'o','MarkerSize',Marker_Size,'color',color_table(5,:),'MarkerFaceColor',color_table(5,:))


X_ticks_2=repmat(0.2,1,length(X_ticks));

glc_bennet=X_ticks;

plot(log2(X_ticks_2),log2(glc_bennet),'o','MarkerSize',Marker_Size,'color',color_table(5,:),'MarkerFaceColor',color_table(5,:))



%% Concentrations biggar experiments

Gal_biggar=[2^-14 0.5 0.2 2];
Glc_biggar=repmat(0.00001,1,length(Gal_biggar));


plot(log2(Gal_biggar),log2(Glc_biggar),'^','MarkerSize',Marker_Size,'color',[0 0 0],'MarkerFaceColor',color_table(2,:));

%%
% Glc_biggar_2=[0 0.075 0.1 0.125 0.15 0.2 0.5];
% Gal_biggar_2=[2 2 2 2 2 2 2];
% 
% 
% plot(log2(Gal_biggar_2),log2(Glc_biggar_2),'^','MarkerSize',Marker_Size,'color',color_table(2,:));

Glc_biggar_3=[0.005 0.02 0.075 0.1 0.125 0.15 0.2 0.5 1]; %;
Gal_biggar_3=repmat(2,1,9);


plot(log2(Gal_biggar_3),log2(Glc_biggar_3),'^','MarkerSize',Marker_Size,'color',[0 0 0],'MarkerFaceColor',color_table(2,:));

%% Concentrations Acar experiments

Glc_acar=[0.1 0.1 0.1 0.1 0.1 0.1];
Gal_acar=[0.001 0.025 0.05 0.1,0.25 0.4];


plot(log2(Gal_acar),log2(Glc_acar),'s','MarkerSize',Marker_Size,'color',[0 0 0],'MarkerFaceColor',color_table(4,:));

%% Concentrations Maier experiments

gal_maier=[0 0.001 0.004 0.01 0.04 1];
glc_maier=repmat(0.02,1,length(gal_maier));


plot(log2(gal_maier),log2(glc_maier),'d','MarkerSize',Marker_Size,'color',[0 0 0],'MarkerFaceColor',color_table(3,:))


%% Concentrations Venturelli experiments


gal_venturelli=[0.002 0.004 0.008 0.02 0.03 0.07 0.15];
glc_venturelli=repmat(0.00001,1,length(gal_venturelli));


plot(log2(gal_venturelli),log2(glc_venturelli),'v','MarkerSize',Marker_Size,'color',[0 0 0],'MarkerFaceColor',color_table(6,:))

%%

xlim([-14.5 2])
ylim([-17.5 1])

Set_fig_RE(hfig2,18,18,22);

%%
filename='Concentrations_Previously_Assayed_5';
export_fig(filename, '-pdf','-transparent','-nocrop');

end
