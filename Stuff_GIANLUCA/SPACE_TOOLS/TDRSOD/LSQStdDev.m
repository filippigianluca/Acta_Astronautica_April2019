function sigma = LSQStdDev(R,n)

sigma = zeros(n,1);

% Covariance
C = LSQCov(n,R);
  
% Standard deviation
for i=1:n;
    sigma(i)=sqrt(C(i,i));
end

