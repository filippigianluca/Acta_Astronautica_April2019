function [f]=quadx_FZ(df_dx,w)

% Quadrature with a fixed nodes and pre-computed integral points
% integration interval

% Federico Zuiani 22/05/2013

% [x,w]=lgwt(N,a,b);
% [df_dx]=fun([x;b]);
 
f=sum(df_dx.*w.');


return