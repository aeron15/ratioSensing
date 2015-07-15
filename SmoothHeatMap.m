function [x,y,s,a,b,a_d,a_u,b_d,b_u] = SmoothHeatMap(M,max,min,cutoff,gal,glu,fit_cuttoff,mid_value)

figure(100)
NM = (M-min)./(max-min);
contour_line = contour((gal),(glu),NM,cutoff);
x = contour_line(1,2:end)';
y = contour_line(2,2:end)';

x = x(find(x~=0));y = y(find(x~=0));
x = x(find(y~=0));y = y(find(y~=0));
x = x(find(y<=fit_cuttoff(1)));y = y(find(y<=fit_cuttoff(1)));
x = x(find(y>=fit_cuttoff(2)));y = y(find(y>=fit_cuttoff(2)));


s = fit(log2(x),log2(y),'x+b','robust','on');
% plot(log2(x),log2(y),'.');hold on;
% plot(log2(x),s(log2(x)));

a =1;% s.a;
b = s(log2(mid_value))/log2(mid_value);
C = confint(s,0.682);
CC = predint(s,log2(mid_value),0.682)
a_d = C(1,1);
a_u = C(2,1);
b_d = CC(1)/log2(mid_value);
b_u = CC(2)/log2(mid_value);

x = contour_line(1,2:end)';
y = contour_line(2,2:end)';

close (gcf)