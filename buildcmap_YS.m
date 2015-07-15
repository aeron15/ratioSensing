function [cmap]=buildcmap_YS(colors)

% Input is a matrix M with size of nX3. each row is a color.

%--------------------------------------------------------------------------


L = size(colors)-1;

ncolors = L(1)

bins=round(255/ncolors);
% diff1=255-bins*ncolors;

vec=zeros(300,3);

vec(1,:) = colors(1,:);


for i=1:ncolors
    beG=(i-1)*bins+1;
    enD=i*bins+1; %beG,enD
    
    vec(beG:enD,1)=linspace(vec(beG,1),colors(i+1,1),bins+1)';
    vec(beG:enD,2)=linspace(vec(beG,2),colors(i+1,2),bins+1)';
    vec(beG:enD,3)=linspace(vec(beG,3),colors(i+1,3),bins+1)';
    
end

cmap=vec(1:bins*ncolors,:);
end %end of buildcmap
