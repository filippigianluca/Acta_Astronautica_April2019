function [f] = f_damper(d,u,par)

C2 = d(1);
T = d(2);
B = u(1);
	

Mu = 0.1;
C1 = 0.1;
%w1 = 100;

Z = sqrt(((B/T)^2*(B^2-1) - (1+Mu)*B^2  -4*C1*C2*B^2/T +1)^2 + 4*( C1*B^3/T^2 + (C2*B^3*(1+Mu) -C2*B)/T -C1*B )^2);

f = (1/Z)*sqrt( (1-(B/T)^2)^2 + 4*(C2*B/T)^2);

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end