function [f] = mv_1_3(d,u,par)
% Multidimensional test function MV3

f1  = sum( d.*u.^2 );
f3  = sum( (5-d).*(1+cos(u)) + (d-1).*(1+sin(u)) );

fs = par.lambda.*[(f1 - 50)/99.8935, (f3 - 13.6569)/2.3431];

f = max(fs);

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
