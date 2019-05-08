function [f] = f_sonc_14(d,u,par)

betas = [.7 .75 1 1.2 ];
sigmas = [.3 .4 1 .4];
peaks0 = [1 1 6 7 8 1 1 6 7 8 1 1 6 7 8 1 1 6 7 8 1 1 6 7 8 ; 
	1 3 8 9.5 2 1 3 8 9.5 2 1 3 8 9.5 2 1 3 8 9.5 2 1 3 8 9.5 2 ; 
	3 1 3 2 5 3 1 3 2 5 3 1 3 2 5 3 1 3 2 5 3 1 3 2 5 ; 
	3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5]; %initially the highest peak
npeaks = 4;

ii = 1:length(d);

umax = d.^ii; % here is where the highest peak will be [0,1]^dim(D)
umax = 5*umax;
uu = 5*u;

for i=1:npeaks
	peak = peaks0(i,:) - [3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5 3 4 1.3 5 5] + umax;
	dist2(i) = sum((uu- peak).^2);
end

f_u = sum(betas.*exp(-dist2./(2.0*sigmas.^2))); %[0,1.2]

f_d = sum(d.^2); % sphere [0,dim(D)]

f = f_d + f_u;


global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end