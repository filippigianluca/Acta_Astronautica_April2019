function [f] = f_rbf_10d(d,u,par)

betas = [.7 .75 1 1.2 1 .6 .5 .2 .4 .1];
sigmas = [.3 .4 1.0 .4 .6 .5 .1 1 .2 .3];
peaks = [1 1 6 7 8 1 1 6 7 8; 1 3 8 9.5 2 1 3 8 9.5 2; 3 1 3 2 5 3 1 3 2 5; 3 4 1.3 5 5 3 4 1.3 5 5; 5 2 9.6 7.3 8.6 5 2 9.6 7.3 8.6; 
	7.5 8 9 3.2 4.6 7.5 8 9 3.2 4.6; 5.7 9.3 2.2 8.4 7.1 5.7 9.3 2.2 8.4 7.1; 5.5 7.2 5.8 2.3 4.5 5.5 7.2 5.8 2.3 4.5 ; 4.7 3.2 5.5 7.1 3.3 4.7 3.2 5.5 7.1 3.3 ; 9.7 8.4 .6 3.2 8.5 9.7 8.4 .6 3.2 8.5];
npeaks = 10;

x = d+u;
dist2 = zeros(1,npeaks);

for i = 1:npeaks
	dist2(i) = sum((x - peaks(i,:)).^2);
end

f = -sum(betas.*exp(-dist2./(2.0*sigmas.^2)));

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end