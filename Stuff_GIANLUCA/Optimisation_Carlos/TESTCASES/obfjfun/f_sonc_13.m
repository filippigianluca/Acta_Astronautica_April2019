function [f] = f_sonc_13(d,u,par)

dx = d(1:2:end);
dy = d(2:2:end);

dx2 = dx.^2;
dy2 = dy.^2;

hd = sum((dx2+dy-11*ones(size(dx))).^2+(dx+dy2-7*ones(size(dx))).^2);

ru = 10*length(u)+sum(u.^2-10*cos(2*pi*u));

f = hd*(1.0-ru/250.0);

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end