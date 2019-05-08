%--------------------------------------------------------------------------
% Initialization with apriori information
% 
% Inputs:
% nEst   Number of estimation parameters
% x      a priori parameters
% P      a priori covariance
% 
% Outputs:
% R, d      
%--------------------------------------------------------------------------
function [] = LSQInit(nEst, x, P)

global R d

% Start the factorization of matrix P. Compute upper triangular
% factor R of P, where P = (R)*(R^T). Proceed backward column
% by column.
for j = nEst:-1:1

% Compute j-th diagonal element.
Sum = 0;
for k = j+2:nEst
    Sum = Sum + R(j,k)*R(j,k);
end
R(j,j) = sqrt(P(j,j)-Sum);

% Complete factorization of j-th column.
for i = j:-1:1
    Sum = 0;
    for k = j+2:nEst
        Sum = Sum + R(i,k)*R(j,k);
    end
    R(i,j) = (P(i,j)-Sum)/R(j,j);
end

end

% Replace R by its inverse R^(-1)
R = InvUpper(R);

% Initialize right hand side
d = R*x;

