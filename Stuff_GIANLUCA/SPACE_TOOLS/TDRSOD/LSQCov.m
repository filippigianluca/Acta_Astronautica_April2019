%--------------------------------------------------------------------------
% 
% Covariance matrix
% 
%--------------------------------------------------------------------------
function T = LSQCov(N,R)

% Calculate the inverse T = R^(-1)  
T = InvUpper(R);

% Replace T by the covariance matrix C=T*T^t
for i=1:N
    for j=i:N
        Sum = 0;
        for k=j:N
            Sum = Sum + T(i,k)*T(j,k);
        end
        T(i,j) = Sum;
        T(j,i) = Sum;
    end
end

