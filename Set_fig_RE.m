function []=Set_fig_RE(fig,axis_fs,title_fs,label_fs)

% Set figure parameters
set(gcf, 'color', [1 1 1]) ;
N = get(fig,'children');


for i = 1:length(N);
set(N(i),'fontsize',axis_fs);

    
%     set(N(i),'fontsize',axis_fs,...
%         'fontname','Times');
    
%     x = get(N(i),'xlabel');set(x,'fontsize',label_fs,'fontname','Cambria')
%     y = get(N(i),'ylabel');set(y,'fontsize',label_fs,'fontname','Cambria')
%     t = get(N(i),'Title');set(t,'fontsize',title_fs,'fontname','Cambria')
     
    x = get(N(i),'xlabel');set(x,'fontsize',label_fs)
    y = get(N(i),'ylabel');set(y,'fontsize',label_fs)
    t = get(N(i),'Title');set(t,'fontsize',title_fs)

end