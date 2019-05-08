%--------------------------------------------------------------------------
%
%  Add a data equation of form Ax=b to the least squares system
%  (performs a row-wise QR transformation using Givens rotations)
% 
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [] = LSQAccumulate(N, A, b, sigma)

global R d

% Weighting
a = A/sigma;  % Normalize A 
b = b/sigma;  % Normalize b

% Construct and apply Givens plane rotation.
for i=1:N
    % Construct the rotation and apply it to
    % eliminate the i-th element of a.
    if ( R(i,i)==0 &&  a(i)==0 )
        c = 1;
        s = 0;
        R(i,i) = 0;
    else
        h = sqrt ( R(i,i)*R(i,i) + a(i)*a(i) );
        if ( R(i,i)<0 )
            h = -h;
        end
        c = R(i,i)/h;
        s = a(i)/h;
        R(i,i) = h;
    end
    
    a(i) = 0;
    
    % Apply the rotation to the remaining elements of a
    for j=i+1:N         
        h       = +c*R(i,j)+s*a(j);
        a(j)    = -s*R(i,j)+c*a(j);
        R(i,j) = h;
    end
    
    % Apply the rotation to the i-th element of d
    h    = +c*d(i)+s*b;
    b    = -s*d(i)+c*b;
    d(i) = h;
end

