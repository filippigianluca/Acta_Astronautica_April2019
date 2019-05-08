function [xi_scaled, c] = chebpolfit(fname,n,a,b,flag_log)

% Input:
%        fname -> function to fit
%        n -> required polynomial degree
% Output:
%        c -> Chebyshev polynomial coefficient (n+1)
%        x -> n+1 Chebyshev nodes

% Marilena Di Carlo, 2015
%% f needs to be computed in correspondence of xi' (that is, xi scaled accorigingly to the interval)
if a == -1 && b == 1
    % Chebyshev nodes - column variable with (n+1) elements
    xi_scaled = cos( (2*(1:n+1)'-1)*pi / (2 * (n+1) ) );
else
    xi_scaled = (a+b)/2 + ((b-a)/2) *cos( (2*(1:n+1)'-1)*pi / (2 * (n+1) ) );
end

for i = 1 : length(xi_scaled)
    % Function values corresponding to Chebyshev nodes
    if flag_log == 1
        y(i,1) = log(feval(fname,xi_scaled(i)));
    else
        y(i,1) = (feval(fname,xi_scaled(i)));
    end
end


%% To compute c, T must be computed as T(xi), not T(xi_scaled)
xi = cos( (2*(1:n+1)'-1)*pi / (2 * (n+1) ) );

% Chebyshev polynomials - 1st column is T0 for all the nodes, 2nd column
% is T1 for all the nodes and so on, for a total of (n+1) columns
% T = [T0 T1 T2 T3 .... Tn]
T = [ones(n+1,1) xi zeros(n+1,n-1)];

% Coefficient (they are n)
c = [sum(y)/(n+1) zeros(1,n)];


for k = 2 : n+1
    
    if k ~= n+1
        
        % Recursive formula for the computation of T
        T(:,k+1) = 2*xi.*T(:,k) - T(:,k-1);
        
    end
    
    c(k) = sum( T(:,k) .* y) * 2 / (n+1);
    

end

