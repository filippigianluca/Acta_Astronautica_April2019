function [f] = f_rbf_2d(d,u,par)

betas = [.7 .75 1 1.2 1];
sigmas = [.3 .4 1.0 .4 .6];
peaks = [1 1 ; 1 3 ; 3 1; 3 4; 5 2];
npeaks = 5;

x = d+u;
dist2 = zeros(1,npeaks);

for i = 1:npeaks
	dist2(i) = sum((x - peaks(i,:)).^2);
end

f = -sum(betas.*exp(-dist2./(2.0*sigmas.^2)));

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end