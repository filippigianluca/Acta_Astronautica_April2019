function [f] = f_soc_14(d,u,par)

ii = 1:length(d);

f = sum(d.^2-((d.^ii)-u).^2);

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end