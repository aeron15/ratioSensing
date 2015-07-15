clear all
close all
clc

syms n(t) s(t) g u N0 S0

z = dsolve(diff(n) == n*g, diff(s) == -n*s*u,n(0)==N0,s(0)==S0);

z.s(1)
t_vec = linspace(0,2,100);
S0_vec = linspace(1,10,10);
for i = 1:length(S0_vec)
S(:,i) = subs(z.s,{u,g,N0,t,S0},{1,1,1,t_vec,S0_vec(i)});
N(:,i) = subs(z.n,{u,g,N0,t,S0},{1,1,1,t_vec,S0_vec(i)});
end
figure(1)
plot(t_vec,N,t_vec,S)

%%
S(30,:)