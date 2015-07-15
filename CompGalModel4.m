% Gal computaintal mdoel
%%
clear all
close all
clc


% Parameters C = [nM] time = [min]

% 1% = 55nM

kf83 = 100;     x(1) = kf83;
kr83 = 150;       x(2) = kr83;
kf84 = 80;     x(3) = kf83;
kr84 = 1300;      x(4) = kr84;

aG1 = 15;       x(5) = aG1;
aG3 = 15;       x(6) = aG3;
aG4 = 2;        x(7) = aG4;
a0G80 = 1*3.1;      x(8) = a0G80;
aG80 = 1*6;    x(9) = aG80;
amig = 1*0.7;       x(10) = amig;

KG1 = 50;        x(11) = KG1;
KG3 = 50;        x(12) = KG3;
KG80 = 50;       x(13) = KG80;

Kmig1 = 50;     x(14) = Kmig1;
Kmig3 = 50;      x(15) = Kmig3;
Kmig4 = 50;     x(16) = Kmig4;

n1 = 3;         x(17) = n1;
n3 = 2;         x(18) = n3;1
n80 = 2;        x(19) = n80;

d = 0.004;
gG1 = d;        x(20) = gG1;
gG3 = d;        x(21) = gG3;
gG4 = d;        x(22) = gG4;
gG80 = d;       x(23) = gG80;
gG33 = d;       x(24) = gG33;
gG83 = d;       x(25) = gG83;
gG84 = d;       x(26) = gG84;
gMig = d;   x(27) = gMig;

% Kgal = 0.02*55;  x(28) = Kgal;% 1mM = 0.02%
Kglu = 150;   x(29) = Kglu;%2^-6% = 0.0156%



kf33 = 70;       x(30) = kf33;
kr33 = 3391;         x(31) = kr33;

vmax_hxt = 200;      x(32) = vmax_hxt;
Khxt_glc = 1;     x(33) = Khxt_glc;
Khxt_gal = 250;     x(34) = Khxt_gal;

vmax_gal2 = 1*1;    x(35) = vmax_gal2;
Kgal2_glc = 5;    x(36) = Kgal2_glc;
Kgal2_gal = 1;    x(37) = Kgal2_gal ;

e=1;

%for r = 1:1;
r=1;
    if r==1
    end
    if r==2
        x(5) = x(5)/2;
    end
    if r==3;
        vmax_gal2 = vmax_gal2/2;
    end
    if r==4
        x(6) = x(6)/2;
    end
    if r==5
        x(7)=x(7)/2;
    end
    if r==6
        x(8) = x(8)/2;
        x(9) = x(9)/2;
    end
    
    %Glu = 0
    load('C:\Users\ys151\Dropbox\gal_paper\Data\s288c\mean_matrix_bfp_yfp');
    aglu_vec = 0;
    agal_vec =  logspace(-3,0.3,50)*55;
    
    gal{1} = [0 2.^[-9:0.5:2]];
    glu{1} = [0 2.^[-7:0.5:0]];
    pcolor(log2(gal{1}),log2(glu{1}),flipud(mean_matrix));
    
    [T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,0,0,e,x),[0 12*60],[0 0 0 0 0 0 0 0 0 0 0]); % Solve ODE
    Y0 = Y(end,:);
    
    for i = 1:length(agal_vec);
        j=1
        agal = agal_vec(i) ;
        aglu = aglu_vec(j);
        
        
        glc_in_hxt = vmax_hxt./(1 + Khxt_glc./aglu + Khxt_glc/Khxt_gal*agal./aglu);
        gal_in_hxt = vmax_hxt./(1 + Khxt_gal./agal + Khxt_gal/Khxt_glc*aglu./agal);
        
        glc_in_gal2 = vmax_gal2./(1 + Kgal2_glc./aglu + Kgal2_glc/Kgal2_gal*agal./aglu);
        gal_in_gal2 = vmax_gal2./(1 + Kgal2_gal./agal + Kgal2_gal/Kgal2_glc*aglu./agal);
        
        gal_in(i) = gal_in_hxt+gal_in_gal2;
        glc_in = glc_in_hxt+glc_in_gal2;
        
%         gal_in(i) = agal;
        [T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,gal_in(i),0,e,x),[0 8*60],Y0); % Solve ODE
        G_vec(i,:) = Y(end,:);
    end
    %
    close all
    
    figure(1)
    plot(log2(agal_vec/55),G_vec(:,[4 5 6 7 1 3 ]),'.-');
    legend({'Gal 4' 'Gal3*' 'Gal 80-3*' 'Gal 80-4' 'Gal 1' 'Gal 3'})%Gal4-80' 'Mig1' 'Gal1' 'Gal80' 'Gal3' 'Gal3-80'})
    Set_fig_YS(figure(1),12,12,18);xlabel('Gal [AU]');ylabel('Response [AU]')
    
    figure(2)
    plot(log2(agal_vec/55),mat2gray(G_vec(:,[1])),'.-');hold on
    plot(log2(gal{1}),mat2gray(mean_matrix(end,:)),'r.-');
    
    figure(3)
    plot(log2(agal_vec/55),gal_in,'.-');hold on;plot(log2(agal_vec/55),agal_vec,'r.-')
    
    %% Double Gradiant
    
    clear gal_in
    
    
    clear f_vec G4_vec RA_vec RE_vec Y YY
    
    gal_res = [-9:0.5:2];
    glc_res = [-7:0.5:0];
    
    agal_vec = [ 2.^gal_res]*55;
    aglu_vec = [ 2.^glc_res]*55;
    
    for i = 1:length(agal_vec);
        for j = 1:length(aglu_vec);
            [i,j]
            agal = agal_vec(i) ;
            aglu = aglu_vec(j);
            
            
            glc_in_hxt = vmax_hxt./(1 + Khxt_glc./aglu + Khxt_glc/Khxt_gal*agal./aglu);
            gal_in_hxt = vmax_hxt./(1 + Khxt_gal./agal + Khxt_gal/Khxt_glc*aglu./agal);
            
            glc_in_gal2 = vmax_gal2./(1 + Kgal2_glc./aglu + Kgal2_glc/Kgal2_gal*agal./aglu);
            gal_in_gal2 = vmax_gal2./(1 + Kgal2_gal./agal + Kgal2_gal/Kgal2_glc*aglu./agal);
            
            gal_in(i,j) = agal; %gal_in_hxt+gal_in_gal2;
            glc_in(i,j) = aglu;%glc_in_hxt+glc_in_gal2;
            
            
            [T Y] = ode15s(@(t,y) CompGalModelEqs(t,y,gal_in(i,j),glc_in(i,j),e,x),[0 8*60],Y0); % Solve ODE
            
            for k = 1:11
                YY(i,j,k) = Y(end,k);
            end
            
            
        end
    end
    
    agal_vec = agal_vec/55;
    aglu_vec = aglu_vec/55;
    
    %%
    close all
    cmap = buildcmap('bwr');
    cmap = 'jet';
    
    figure(1)
    subplot(1,2,1)
    surf(log2(agal_vec),log2(aglu_vec),glc_in');hold on;view(0,90);
    axis('square');
    xlabel('Gal');ylabel('Glu');title('GLC in');
    colormap(cmap);
    subplot(1,2,2)
    surf(log2(agal_vec),log2(aglu_vec),gal_in');hold on;view(0,90);
    axis('square');
    xlabel('Gal');ylabel('Glu');title('GAL in');
    colormap(cmap);
    
    figure(3)
    surf(log2(agal_vec),log2(aglu_vec),YY(:,:,4)');hold on;view(0,90);
    axis('square');
    xlabel('Gal');ylabel('Glu');title('Gal 4');
    colormap(cmap);
    
    figure(4)
    surf(log2(agal_vec),log2(aglu_vec),(YY(:,:,1)'));hold on;view(0,90);
    axis('square');
    xlabel('Gal');ylabel('Glu');title('Gal 1')
    colormap(cmap);
    xlim([-9 1])
    
    figure(5)
    surf(log2(agal_vec),log2(aglu_vec),YY(:,:,11)');hold on;view(0,90);
    axis('square');
    xlabel('Gal');ylabel('Glu');title('Mig1')
    colormap(cmap);
    
    figure(6)
    pcolor(log2(gal{1}),log2(glu{1}),flipud(mean_matrix));
    axis('square');
    xlabel('Gal');ylabel('Glu');title('GAl1 experiment')
    colormap(cmap);
    % Compute model slope
    
    fit_cuttoff = [2^-1 2^-6];
    mid_value = 2^-4;
    cutoff=0.2;
    i=1
    [xx,y,s,a(i),b(i),a_d(i),a_u(i),b_d(i),b_u(i)] = SmoothHeatMap(YY(:,:,1)',max(max(YY(:,:,1)')),min(min(YY(:,:,1)')),cutoff,agal_vec,aglu_vec,fit_cuttoff,mid_value);
    
    figure(7);
    plot(log2(xx),log2(y),'ok','markersize',1.5);hold on;
    plot(log2(xx),s(log2(xx)),'k');hold on;
    
    slope_front(r) = a
    cutoff_front(r) = b
    xxx{r} = xx;yy{r}=y;ss{r} = s(log2(xx));
    
%end
%%
figure(8)
subplot(2,1,1)
bar(slope_front)
subplot(2,1,2)
bar(cutoff_front)


figure(9)
for i = 1:6
    plot(log2(xxx{i}),log2(yy{i}),'o','markersize',1.5);hold all;
end

for i = 1:6
    
    plot(log2(xxx{i}),ss{i});hold all;
end

legend([ {'WT' 'GAl1' 'Gal2' 'Gal3' 'Gal4' 'Gal80'}])

figure(10)
subplot(2,1,1);
plot(slope_front(2:6)/slope_front(1),'o');hold on;refline(0,1)
subplot(2,1,2)
plot(slope_front(2:6)/slope_front(1),'o');hold on;refline(0,1)
