function [c,ceq] = DirTranscFunctionFMINCONGLConst(y,inCond,fCond,N,h,mu,uMax)
%
% Function DirTranscFunctionFMINCONGLConst: constraints of the direct transcription of the low thrust
% optimal control problem. Implements a 5th order Gauss-Lobatto scheme
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
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the constraints

% Extract the state and control
r = y(1:3,1:2*N-1);
v = y(4:6,1:2*N-1);
u = y(7:9,1:2*N-1);

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

% Transcription constraints (5th order Gauss-Lobatto scheme)
for i=1:2:2*N-3
    % Calculate rdotk and vdotk, the derivativves at time k (i)
    rdotK = v(1:3,i);
    vdotK = -mu*r(1:3,i)/norm(r(1:3,i)).^3 + u(1:3,i);
    
    %%%%%%%%%%%%%% Define ceqR and ceqV with 5th order Gauss-Lobatto scheme %%%%%%%%%%%%%%%%%%%
    
    %Calculate rdot and vdot at k+2 (because of the notation selected
    %correspondent to an higher discretized interval)
    rdotK1 = v(1:3,i+2);
    vdotK1 = -mu*r(1:3,i+2)/norm(r(1:3,i+2)).^3 + u(1:3,i+2);
    
    %Calculate rdot and vdot at C, middle point between k and k+1
    %(available in the variable vector due to the increased
    %discretisation)
    rdotC = v(1:3,i+1);
    vdotC = -mu*r(1:3,i+1)/norm(r(1:3,i+1)).^3 + u(1:3,i+1);
    
    % Find the values of the collocation points x1,x2
    r1 = (1/686)*( (39*sqrt(21) + 231)*r(1:3,i) + 224*r(1:3,i+1) +...
        + (-39*sqrt(21) + 231)*r(1:3,i+2) + h*( (3*sqrt(21) + 21)*rdotK -16*sqrt(21)*rdotC +...
        + (3*sqrt(21) - 21)*rdotK1 ) );
    
    v1 = (1/686)*( (39*sqrt(21) + 231)*v(1:3,i) + 224*v(1:3,i+1) +...
        + (-39*sqrt(21) + 231)*v(1:3,i+2) + h*( (3*sqrt(21) + 21)*vdotK -16*sqrt(21)*vdotC +...
        + (3*sqrt(21) - 21)*vdotK1 ) );
    
    r2 = (1/686)*( (-39*sqrt(21) + 231)*r(1:3,i) + 224*r(1:3,i+1) +...
        + (39*sqrt(21) + 231)*r(1:3,i+2) + h*( (-3*sqrt(21) + 21)*rdotK +16*sqrt(21)*rdotC +...
        + (-3*sqrt(21) - 21)*rdotK1 ) );
    
    v2 = (1/686)*( (-39*sqrt(21) + 231)*v(1:3,i) + 224*v(1:3,i+1) +...
        + (39*sqrt(21) + 231)*v(1:3,i+2) + h*( (-3*sqrt(21) + 21)*vdotK +16*sqrt(21)*vdotC +...
        + (-3*sqrt(21) - 21)*vdotK1 ) );
    
    % Estimate the control at the collocation points with a linear interpolation
    % (no other information is given)
    du1 = u(1:3,i+1) - u(1:3,i);
    u1 = u(1:3,i+1) - sqrt(3/7)*du1;
    
    du2 = u(1:3,i+2) - u(1:3,i+1);
    u2 = u(1:3,i+1) + sqrt(3/7)*du2;
    
    %Calculate rdot and vdot at collocation points x1,s2
    rdot1 = v1;
    vdot1 = -mu*r1/norm(r1).^3 + u1;
    
    rdot2 = v2;
    vdot2 = -mu*r2/norm(r2).^3 + u2;
    
    %Write constraints
    ceqR1(3*(i-1)+1:3*(i-1)+3) = (1/360)*( (32*sqrt(21) + 180)*r(1:3,i) - 64*sqrt(21)*r(1:3,i+1) +...
        + (32*sqrt(21) - 180)*r(1:3,i+2) + h*( (9 + sqrt(21))*rdotK + 98*rdot1 + 64*rdotC +...
        + (9 - sqrt(21))*rdotK1 ) );
    
    ceqV1(3*(i-1)+1:3*(i-1)+3) = (1/360)*( (32*sqrt(21) + 180)*v(1:3,i) - 64*sqrt(21)*v(1:3,i+1) +...
        + (32*sqrt(21) - 180)*v(1:3,i+2) + h*( (9 + sqrt(21))*vdotK + 98*vdot1 + 64*vdotC +...
        + (9 - sqrt(21))*vdotK1 ) );
    
    ceqR2(3*(i-1)+1:3*(i-1)+3) = (1/360)*( (-32*sqrt(21) + 180)*r(1:3,i) + 64*sqrt(21)*r(1:3,i+1) +...
        + (-32*sqrt(21) - 180)*r(1:3,i+2) + h*( (9 - sqrt(21))*rdotK + 98*rdot2 + 64*rdotC +...
        + (9 + sqrt(21))*rdotK1 ) );
    
    ceqV2(3*(i-1)+1:3*(i-1)+3) = (1/360)*( (-32*sqrt(21) + 180)*v(1:3,i) + 64*sqrt(21)*v(1:3,i+1) +...
        + (-32*sqrt(21) - 180)*v(1:3,i+2) + h*( (9 - sqrt(21))*vdotK + 98*vdot2 + 64*vdotC +...
        + (9 + sqrt(21))*vdotK1 ) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------------------------------------------------------------------
    % Add constraints on the modulus of control acceleration
    %--------------------------------------------------------------------------
    % c(i) = norm(u(:,i)) - uMax;
end

% c(N) = norm(u(:,N)) - uMax;
ceq = [ceq ceqR1 ceqV1 ceqR2 ceqV2];

end