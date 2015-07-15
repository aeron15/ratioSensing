%% Parse heat mapm for publicatoin and data analyissi


function [D,M] = ParseHeatmapMat(E)

x = padarray(E,[2 2],nan);
x(2,:) = x(3,:);x(3,:)=nan;
x(:,2) = x(:,3);x(:,3)=nan;
x(2,2)=nan;
M = x(2:end-1,2:end-1);
D = x(4:end-2,4:end-2);

