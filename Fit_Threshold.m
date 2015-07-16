function [x,y,s] = Fit_Threshold(M,max,min,cutoff,gal,glu,fit_cuttoff,mid_value)
%function [x,y,s,a,b,a_d,a_u,b_d,b_u] = Fit_Threshold(M,max,min,cutoff,gal,glu,fit_cuttoff,mid_value)

figure(100)
NM = (M-min)./(max-min);
contour_line = contour((gal),(glu),NM,cutoff);
x = contour_line(1,2:end)';
y = contour_line(2,2:end)';

x = x(find(x~=0));y = y(find(x~=0));
x = x(find(y~=0));y = y(find(y~=0));
x1=x;
y1=y;

x = x(find(y<=fit_cuttoff(1)));y = y(find(y<=fit_cuttoff(1)));
x = x(find(y>=fit_cuttoff(2)));y = y(find(y>=fit_cuttoff(2)));


s = fit(log2(y),log2(x),'a*x+b','robust','on');

close (gcf)