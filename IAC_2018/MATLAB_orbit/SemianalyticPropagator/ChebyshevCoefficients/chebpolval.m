function P = chebpolval(x, c, n, a, b)

% Marilena Di Carlo, 2015


if a == -1 && b == 1 
else
    x = -1 + 2* ((x-a)/(b-a));
end

% Compute T for input x
% Chebyshev polynomials - 1st column is T0 for all the nodes, 2nd column
% is T1 for all the nodes and so on, for a total of (n+1) columns
% T = [T0 T1 T2 T3 .... Tn]
T = [ones(length(x),1) x zeros(length(x), n-1)];


for k = 2 : n+1
    
    if k ~= n+1
        
        % Recursive formula for the computation of T
        T(:,k+1) = 2*x.*T(:,k) - T(:,k-1);
        
    end
   

end


P = c(1) ;

for k = 2 : n+1
    P = P + c(k) * T(:,k);    
end


end