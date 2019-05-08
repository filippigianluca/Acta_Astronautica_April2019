function [P] = nfe2pos(n_FE, n_int)

prodotto = prod(n_int);  
dim=length(n_int);

A =1;
for j = dim:-1:2
    A = A*n_int(j);
    P(j) = ceil(n_FE*A/prodotto);
    n_FE = n_FE-prodotto/A*(P(j)-1);
end
P(1) = n_FE;

end