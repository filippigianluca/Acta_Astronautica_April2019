%--------------------------------------------------------------------------
% Solve the LSQ problem for vector x[] by backsubstitution
%--------------------------------------------------------------------------
function x = LSQSolve(N, R, d)
  
% Check for singular matrix 
for i=1:N
    if ( R(i,i) == 0 )
        error(' ERROR: Singular matrix R in LSQSolve()');
    end
end

% Solve Rx=d for x_n,...,x_1 by backsubstitution
x(N) = d(N) / R(N,N);
for i=N-1:-1:1
    Sum = 0;
    for j=i+1:N
        Sum = Sum + R(i,j)*x(j);
    end
    x(i) = ( d(i) - Sum ) / R(i,i);
end

