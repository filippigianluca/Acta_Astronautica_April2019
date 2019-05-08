function [X,Y,M] = porkchop_construction(name)

A = importdata(name);

[X,i] = unique(A.data(:,1));
% X = X(1:end-1);
l_dep = length(X);

[Y,j] = unique(A.data(:,2));
l_TOF = length(Y);

nr = max(A.data(:,3));

DV = A.data(:,4) + A.data(:,5).*0.04;
M = max(DV).*ones(l_dep,l_TOF);

for ii = 1:l_dep
    for jj = 1:l_TOF
        
        M(ii,jj) = min(DV( i(ii) + j(jj) - 1 : i(ii) + j(jj) + nr - 2 ));
        
    end
end

M = M';

end
