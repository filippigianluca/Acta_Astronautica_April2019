function [f] = f_soc_13(d,u,par)

ux = u(1:2:end);
uy = u(2:2:end);

ux2 = ux.^2;
uy2 = uy.^2;

h = sum((ux2+uy-11*ones(size(ux))).^2+(ux+uy2-7*ones(size(ux))).^2)/(1e3*length(ux));

fd = 10*length(d)+sum(d.^2-10*cos(2*pi*d));

f = fd*(1.0-h); 

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end