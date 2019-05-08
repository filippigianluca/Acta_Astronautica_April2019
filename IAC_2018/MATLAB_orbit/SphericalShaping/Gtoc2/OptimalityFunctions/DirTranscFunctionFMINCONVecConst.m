function [ceq] = DirTranscFunctionFMINCONVecConst(y,inCond,fCond,N,h,mu,uMax)
%
% Function DirTranscFunctionFMINCONVecConst: constraints of the direct transcription of the low thrust
% optimal control problem
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
% Note: written for use with OPTI, so the constraints have to be written in
% a vectorial way instead than matricial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the constraints

% Extract the state and control
r = y(1:3*N);
v = y(3*N+1:6*N);
u = y(6*N+1:9*N);

% Write down the initial and final conditions
r0 = inCond.r0;
v0 = inCond.v0;
rf = fCond.rf;
vf = fCond.vf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear inequalities
%c = [];

% Nonlinear equalities

% Initial and final constraints on the state
ceq(1:3) = r(1:3) - r0';
ceq(4:6) = r(3*(N-1)+1:3*N) - rf';
ceq(7:9) = v(1:3) - v0';
ceq(10:12) = v(3*(N-1)+1:3*N) - vf';

% Transcription constraints
for i=1:N-1
    % Calculate rdot and vdot
    rdotK = v(3*(i-1)+1:3*(i-1)+3);
    vdotK = -mu*r(3*(i-1)+1:3*(i-1)+3)/norm(r(3*(i-1)+1:3*(i-1)+3)).^3 + u(3*(i-1)+1:3*(i-1)+3);
    
%%%%%%%%%%%%%% Define ceqR and ceqV with Euler scheme %%%%%%%%%%%%%%%%%%%%%

%       ceqR(3*(i-1)+1:3*(i-1)+3) = r(3*i+1:3*i+3) - r(3*(i-1)+1:3*(i-1)+3) - h*rdotK;
%       ceqV(3*(i-1)+1:3*(i-1)+3) = v(3*i+1:3*i+3) - v(3*(i-1)+1:3*(i-1)+3) - h*vdotK;
    
%%%%%%%%%%%%%% Define ceqR and ceqV with Hermite scheme %%%%%%%%%%%%%%%%%%% 
   
 % Calculate rdot and vdot at k+1
    rdotK1 = v(3*i+1:3*i+3);
    vdotK1 = -mu*r(3*i+1:3*i+3)/norm(r(3*i+1:3*i+3)).^3 + u(3*i+1:3*i+3);

%Find the intermediate values of r,v,u (needed only for Hermite
%interpolation case)

    rc = (r(3*(i-1)+1:3*(i-1)+3) + r(3*i+1:3*i+3))/2 + (h/8)*(rdotK - rdotK1);
    vc = (v(3*(i-1)+1:3*(i-1)+3) + v(3*i+1:3*i+3))/2 + (h/8)*(vdotK - vdotK1);
    uc = (u(3*(i-1)+1:3*(i-1)+3) + u(3*i+1:3*i+3))/2;


 % Calculate rdot and vdot at c = (k + k+1)/2
    rdotC = vc;
    vdotC = -mu*rc/norm(rc).^3 + uc;

 % Write constraints
    ceqR(3*(i-1)+1:3*(i-1)+3) = r(3*(i-1)+1:3*(i-1)+3) - r(3*i+1:3*i+3) + (h/6)*(rdotK + 4*rdotC + rdotK1);
    ceqV(3*(i-1)+1:3*(i-1)+3) = v(3*(i-1)+1:3*(i-1)+3) - v(3*i+1:3*i+3) + (h/6)*(vdotK + 4*vdotC + vdotK1);
    
%%%%%%%%%%%%%% Define ceqR and ceqV with Crank Nicholson scheme %%%%%%%%%%%

%     % Calculate rdot and vdot at k+1
%     rdotK1 = v(3*i+1:3*i+3);
%     vdotK1 = -mu*r(3*i+1:3*i+3)/norm(r(3*i+1:3*i+3)).^3 + u(3*i+1:3*i+3);

%     % Write constraints
%     ceqR(3*(i-1)+1:3*(i-1)+3) = r(3*(i-1)+1:3*(i-1)+3) - r(3*i+1:3*i+3) - (h/2)*(rdotK + rdotK1);
%     ceqV(3*(i-1)+1:3*(i-1)+3) = v(3*(i-1)+1:3*(i-1)+3) - v(3*i+1:3*i+3) - (h/2)*(vdotK + vdotK1);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Add constraints on the modulus of control acceleration
%--------------------------------------------------------------------------
% c(i) = norm(u(3*(i-1)+1:3*(i-1)+3)) - uMax;
end

% c(N) = norm(u(3*(N-1)+1:3*N)) - uMax;

ceq = [ceq ceqR ceqV];

end