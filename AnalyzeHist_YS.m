function [prec_const,prec_adaptive,prec_area] = AnalyzeHist_YS(y,x,th_const,off_peak,off_dist_y,off_dist_x);



p = y./sum(y);

if(length(find(x>=th_const)))>1
    ind= find(x>=th_const);
    prec_const = trapz(x(ind),y(ind))/trapz(x,y);
    prec_adaptive=prec_const;
else
    prec_const = 0;
    prec_adaptive=0;
end

% set up a common x 

max_x = max([max(x),max(off_dist_x)]);
min_x = min([min(x),min(off_dist_x)]);
xx = linspace(min_x,max_x,150);
% 
yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
yy(isnan(yy))=0;
yy_off = interp1(off_dist_x,off_dist_y,xx);
yy_off(yy_off<0)=0;
yy_off(isnan(yy_off))=0;
% xx = off_dist_x;

pp = yy/trapz(xx,yy);
p_off =yy_off/trapz(xx,yy_off);
diff_p = pp-p_off;

[val,ind] = max(p_off);
% xx(ind)
diff_p(diff_p<0)=0;
diff_p(1:ind) = 0;

% plot(xx,p_off);hold on;
% plot(xx,diff_p,'r')
prec_area = trapz(xx,diff_p);
% plot(xx,p_off,'r',xx,pp,'k',xx,diff_p)