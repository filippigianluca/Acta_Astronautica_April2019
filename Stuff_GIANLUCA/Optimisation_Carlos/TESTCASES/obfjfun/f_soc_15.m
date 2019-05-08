function [f] = f_soc_15(d,u,par)

f0 = 0.5*sum(d.^4-16*d.^2+5*d)/39.16599;

umax = sin(d); %[-1,1]

udist = u-umax; %[-2,2]

correction = sum(cos(udist) - ones(size(u))); %[-1.41,0] and convex

f = f0+correction;

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end