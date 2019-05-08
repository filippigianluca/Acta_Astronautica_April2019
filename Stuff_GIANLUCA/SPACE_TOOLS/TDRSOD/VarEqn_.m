%--------------------------------------------------------------------------
%
% VarEqn_: Computes the variational equations, i.e. the derivative of the
%          state vector and the state transition matrix
%
% Inputs:
%   t           Time since epoch AuxParam.Mjd_0 in [s]
%   yPhiS       (6+36+6)-dim vector comprising the state vector (y), the
%               state transition matrix (Phi) and the sensitivity matrix
%               in column wise storage order
%
% Output:
%   yPhiSp      Derivative of yPhiS
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function yPhiSp = VarEqn_ (t, yPhiS)

% State vector components
r = yPhiS(1:3);
v = yPhiS(4:6);
Phi = zeros(6);
S = zeros(6,2);
dfdy = zeros(6);
dfdp = zeros(6,2);
yPhiSp = zeros(size(yPhiS));

% State transition matrix
for j=1:6
    Phi(:,j) = yPhiS(6*j+1:6*j+6);
end

% Sensitivity matrix 
for j=7:8
    S(:,j-6) = yPhiS(6*j+1:6*j+6);
end

% Acceleration and gradient
[a, G, dadCD, dadCR] = Accel_Grad(t, [r;v]);

% Time derivative of state transition matrix
for i=1:3
    for j=1:3
        dfdy(i  ,j) = 0;                    % dv/dr(i,j)
        dfdy(i+3,j) = G(i,j);               % da/dr(i,j)        
        if (i==j)
            dfdy(i  ,j+3) = 1; % dv/dv(i,j)
        else
            dfdy(i  ,j+3) = 0; % dv/dv(i,j)
        end               
        dfdy(i+3,j+3) = 0;                  % da/dv(i,j)
    end
end

Phip = dfdy*Phi;

% Time derivative of sensitivity matrix
for i=1:3
    dfdp(i  ,1) = 0;                      % dv/dCD(i)
    dfdp(i+3,1) = dadCD(i);               % da/dCD(i)
    dfdp(i  ,2) = 0;                      % dv/dCR(i)
    dfdp(i+3,2) = dadCR(i);               % da/dCR(i)
end

Sp = dfdy*S + dfdp;

% Derivative of combined state vector and state transition matrix
for i=1:3
    yPhiSp(i)   = v(i);                    % dr/dt(i)
    yPhiSp(i+3) = a(i);                    % dv/dt(i)
end

for i=1:6
    for j=1:6
        yPhiSp(6*(j)+i  ) = Phip(i,j);     % dPhi/dt(i,j)
    end
end

for i=1:6
    for j=7:8
        yPhiSp(6*(j)+i  ) = Sp(i,j-6);     % dS/dt(i,j)
    end
end

