function predict_responses(gal_final,glc_final,gal_m,glc_m,D,M)

%Predict Responses by multiplying the single input responses

%% single response

% gal_r_area = E_area{1}(1,:);
% glc_r_area= E_area{1}(:,end-2); % 2% GAL

cmap=cbrewer('seq', 'YlOrRd', 25);


gal_final=gal_final(2:end);
glc_final=glc_final(2:end);
%%
gal_r_area = D{1}(1,:);
glc_r_area= D{1}(:,end-2); % 2% GAL

% figure;
% plot(glc_r_area)
% figure;
% plot(gal_r_area)
%%
fig4=figure(4);
set(fig4,'Position',[ 1   179   495   627])

subplot(2,1,1)
plot(log2(gal_final),gal_r_area,'o-k','linewidth',2);hold on;
xlim([-9.5 2.5]);set(gca,'xtick',[-9:3:2,2]);
ylim([0 1])
%xlabel('Galactose [log_{2}(%)]')
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 .02 0])

subplot(2,1,2)
plot(log2(glc_final),glc_r_area,'o-k','linewidth',2);hold on;
%xlim([-9.5 0.5]);set(gca,'xtick',[-9:3:0]);
xlim([-10 1.5]);set(gca,'xtick',[-10:3:1,1]);
ylim([0 1])

%xlabel('Glucose [log_{2}(%)]')
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 .02 0])

Set_fig_RE(figure(4),16,16,20)

% filename=['output/Single_responses_WT'];
% export_fig(filename, '-eps','-transparent','-nocrop');

%%
for i = 1:length(gal_r_area)
    for j = 1:length(glc_r_area)
        R(j,i) = gal_r_area(i)*glc_r_area(j);
    end
end
[D_R,M_R] = ParseHeatmapMat(R);
%%

fig1=figure(1)
scrsz = get(0,'ScreenSize');
set(fig1,'Position',[1 scrsz(4)*0.45 scrsz(3)*0.45 scrsz(4)])

h = pcolor(M_R);
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-9:1:0)])
%xlabel('Galactose [log_{2}(%)]')
%ylabel('Glucose [log_{2}(%)]')
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 0.5 0])
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0 0.5 0])

%hc=colorbar;
%ylabel(hc,'fraction of induced cells')

colormap(cmap);
Set_fig_RE(figure(1),18,18,18)

% filename=['Predicted_response_WT'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
%%
fig2=figure(2)
scrsz = get(0,'ScreenSize');
set(fig2,'Position',[1 scrsz(4)*0.45 scrsz(3)*0.45 scrsz(4)])

h = pcolor(M{1});
set(h,'edgecolor','none');
axis square
set(gca,'xtick',[1  3 5:2:25]+0.5,'xticklabel',['No', mat2cell(-9:1:2)]);
set(gca,'ytick',[1  3 5:2:22]+0.5,'yticklabel',['No', mat2cell(-9:1:0)])
% xlabel('Galactose [log_{2}(%)]')
% ylabel('Glucose [log_{2}(%)]')

xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 0.2 0])
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0 0.2 0])

colormap(cmap)
Set_fig_RE(figure(2),18,18,18)

% filename=['Measured_response_WT'];
% export_fig(filename, '-pdf','-transparent','-nocrop');

end
