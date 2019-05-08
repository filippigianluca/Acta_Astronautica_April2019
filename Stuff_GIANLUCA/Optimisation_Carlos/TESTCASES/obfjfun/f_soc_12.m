function [f] = f_soc_12(d,u,par)

d1 = d(1:2:end);
d2 = d(2:2:end);
d12 = d1.^2;
d22 = d2.^2;
u1 = u(1:2:end);
u2 = u(2:2:end);

f = sum(d1.^4 + 2*d22 + (d22).*sin(d22) -(u1-sqrt(d12+d22)).^2 - u2.*(3*ones(size(u1))+d1-2*d2));

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end