%--------------------------------------------------------------------------
%
% InvUpper: Inversion of an upper right triangular matrix
%
% Input:
%   R    Upper triangular square matrix
%
% Output:
%   T    Inverse of R 
%
% Note:
%   This function may be called with the same actual parameter for R and T
%
%--------------------------------------------------------------------------
function [T] = InvUpper(R)

[N, M] = size(R);   % Dimension of R

% Check diagonal elements

if (M ~= N)
    error(' ERROR: Incompatible shapes in InvUpper\n');
    exit(1);
end

for i=1:N
    if ( R(i,i) == 0 )
        error(' ERROR: Singular matrix in InvUpper\n');
        exit(1);
    else
        % Compute the inverse of i-th diagonal element.
        T(i,i) = 1/R(i,i);
    end
end

% Calculate the inverse T = R^(-1)
for i=1:N-1 
    for j=i+1:N
        Sum = 0;
        for k=i:j-1
            Sum = Sum + T(i,k)*R(k,j);
        end
        T(i,j) = -T(j,j)*Sum;
    end
end

