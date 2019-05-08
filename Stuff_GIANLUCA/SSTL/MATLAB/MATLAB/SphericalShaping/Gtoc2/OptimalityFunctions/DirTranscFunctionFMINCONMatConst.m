function [c,ceq] = DirTranscFunctionFMINCONMatConst(y,inCond,fCond,N,h,mu,uMax)
%
% Function DirTranscFunctionFMINCONMatConst: constraints of the direct transcription
% of the low thrust optimal control problem. Different schemes are
% implemented for the constraints (Euler, Hermite, Runge Kutta, Crank
% Nicholson)
% Standard form:
% c = ...     % Compute nonlinear inequalities at y.
% ceq = ...   % Compute nonlinear equalities at y.
% Vector y made of state and control
%
% INPUT
% y: array comprising the state and the control at each moment of time
% y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
% a 3x1 vector (velocity) and u is a 3x1 vector (control)
% It is a matrix [3x3n x 1]
% inCond,fCond: struct with initial and final conditions
% N: number of node points
% h: time step
% mu: gravitational parameter
%
% OUTPUT
% c: nonlinear inequalities. There are no inequalities in this case.
% ceq: nonlinear equalities. The inequalities are 2 + (N-1). The first two
% are the initial and final state constraints. The other (N-1) are the
% transcription constraints. It is an array of size [3x3x2 + 3x3x(N-1)x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the constraints

% Extract the state and control
r = y(1:3,1:N);
v = y(4:6,1:N);
u = y(7:9,1:N);

% Write down the initial and final conditions
r0 = inCond.r0;
v0 = inCond.v0;
rf = fCond.rf;
vf = fCond.vf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear inequalities
% AU2km = 149597870.66;
% T_sid = 86164.1004;  
% Umax = 1.4*1e-7*(T_sid^2/AU2km);

c = [];

% Nonlinear equalities

% Initial and final constraints on the state
ceq(1:3) = r(1:3,1)' - r0;
ceq(4:6) = r(1:3,end)' - rf;
ceq(7:9) = v(1:3,1)' - v0;
ceq(10:12) = v(1:3,end)' - vf;

% Transcription constraints
for i=1:N-1
    % Calculate rdotk and vdotk, the derivativves at time k (i)
    rdotK = v(1:3,i);
    vdotK = -mu*r(1:3,i)/norm(r(1:3,i)).^3 + u(1:3,i);
    
%%%%%%%%%%%%%% Define ceqR and ceqV with Euler scheme %%%%%%%%%%%%%%%%%%%%%

%       ceqR(3*(i-1)+1:3*(i-1)+3) = r(1:3,i+1) - r(1:3,i) - h*rdotK;
%       ceqV(3*(i-1)+1:3*(i-1)+3) = v(1:3,i+1) - v(1:3,i) - h*vdotK;
    
%%%%%%%%%%%%%% Define ceqR and ceqV with Hermite scheme %%%%%%%%%%%%%%%%%%% 

     %Calculate rdot and vdot at k+1
        rdotK1 = v(1:3,i+1);
        vdotK1 = -mu*r(1:3,i+1)/norm(r(1:3,i+1)).^3 + u(1:3,i+1);
    
    %Find the intermediate values of r,v,u (needed only for Hermite
    %interpolation case)
        rc = (r(1:3,i) + r(1:3,i+1))/2 + (h/8)*(rdotK - rdotK1);
        vc = (v(1:3,i) + v(1:3,i+1))/2 + (h/8)*(vdotK - vdotK1);
        uc = (u(1:3,i) + u(1:3,i+1))/2;
    
     %Calculate rdot and vdot at c = (k + k+1)/2
        rdotC = vc;
        vdotC = -mu*rc/norm(rc).^3 + uc;
    
     %Write constraints
        ceqR(3*(i-1)+1:3*(i-1)+3) = r(1:3,i) - r(1:3,i+1) + (h/6)*(rdotK + 4*rdotC + rdotK1);
        ceqV(3*(i-1)+1:3*(i-1)+3) = v(1:3,i) - v(1:3,i+1) + (h/6)*(vdotK + 4*vdotC + vdotK1);

%%%%%%%%%%%%%% Define ceqR and ceqV with Runge-Kutta scheme %%%%%%%%%%%%%%%%%%% 

%     %Find the intermediate values of u (needed only for Runge-Kutta
%     %interpolation case)
%         uc = (u(1:3,i) + u(1:3,i+1))/2;
% 
%      % Calculate k1
%         k1_r = h*rdotK; 
%         k1_v = h*vdotK; 
% 
%      % Calculate k2
%         k2_r = h*(rdotK + k1_v/2);
%         k2_v = h*(-mu*(r(1:3,i) + k1_r/2)/norm(r(1:3,i) + k1_r/2 ).^3 + uc);
% 
%      % Calculate k3
%         k3_r = h*(rdotK + k2_v/2);
%         k3_v = h*(-mu*(r(1:3,i) + k2_r/2)/norm(r(1:3,i) + k2_r/2).^3 + uc);
% 
%      % Calculate k4
%         k4_r = h*(rdotK + k3_v);
%         k4_v = h*(-mu*(r(1:3,i) + k3_r)/norm(r(1:3,i) + k3_r).^3 + u(1:3,i+1));
%     
%      %Write constraints
%         ceqR(3*(i-1)+1:3*(i-1)+3) = r(1:3,i) - r(1:3,i+1) + (1/6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
%         ceqV(3*(i-1)+1:3*(i-1)+3) = v(1:3,i) - v(1:3,i+1) + (1/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);

%%%%%%%%%%%%%% Define ceqR and ceqV with Crank Nicholson scheme %%%%%%%%%%%

%     % Calculate rdot and vdot at k+1
%     rdotK1 = v(1:3,i+1);
%     vdotK1 = -mu*r(1:3,i+1)/norm(r(1:3,i+1)).^3 + u(1:3,i+1);

%     % Write constraints
%     ceqR(3*(i-1)+1:3*(i-1)+3) = r(1:3,i+1) - r(1:3,i) - (h/2)*(rdotK + rdotK1);
%     ceqV(3*(i-1)+1:3*(i-1)+3) = v(1:3,i+1) - v(1:3,i) - (h/2)*(vdotK + vdotK1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Add constraints on the modulus of control acceleration
%--------------------------------------------------------------------------
% c(i) = norm(u(:,i)) - uMax;
end

% c(N) = norm(u(:,N)) - uMax;
ceq = [ceq ceqR ceqV];

end