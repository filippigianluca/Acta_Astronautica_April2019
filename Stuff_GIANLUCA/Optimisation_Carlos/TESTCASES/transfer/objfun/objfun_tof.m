function [tof] = obfjun_tof(d,u,par)
global nfevalglobal;
nfevalglobal = nfevalglobal + 1;
	tof = d(2)+u(2);
	
return