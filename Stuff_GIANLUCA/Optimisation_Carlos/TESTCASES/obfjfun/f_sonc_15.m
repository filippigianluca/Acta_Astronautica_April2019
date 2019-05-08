function [f] = f_sonc_15(d,u,par)

f0 = 0.5*sum(d.^4-16*d.^2+5*d)/39.16599; % styblinski-tang function in d

umax = 5*sin(d); %[-10,10]

dist = sqrt(sum((5*u-umax).^2));

fu = cos(2*dist)/(10+2*dist)-0.1;

f = f0+fu;

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end