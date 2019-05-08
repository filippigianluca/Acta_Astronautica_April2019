function [lam,mu] = cubesat_failure_rate(d,u, lam0, mu0)

if nargin<3
    lam0 = 1;
end

if nargin<4
    mu0 = 3;
end

ad = [0
    0.0033
    -0.11667
    0.1
    -0.1
    0.2
    0
    0
    0
    %0
    0.075
    0.2
    -0.4];

au = [13.33333
    0
    0
    0
    0
    -0.0175
    0
    0
    0.33333
    0.2
    0.75
    0
%     0.00667
    0.015
    0.00333
    0.003
    0
    0
    0.04
    0
    0.01];

bd = [1
    0.8
    2.01667
    0.95
    1.05
    0.9
    1
    1
    1
    %1
    0.825
    0.9
    1.2];

bu = [0.833333
    1
    1
    1
    1
    1.025
    1
    1
    0.816667
    0.8
    0.775
    1  
%     1
    1
    1
    0.85
    1
    1
    0.8
    1  
    0.95];

dlam1 = ad.*d'+bd;
dlam2 = au.*u'+bu;

%safety trick  - not neccessary if we stay in designed margins for d and u
dlam1 = abs(dlam1);
dlam2 = abs(dlam2);

lam = lam0*prod(dlam1)*prod(dlam2);
mu = mu0;

end
