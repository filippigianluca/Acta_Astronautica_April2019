function [r,v,u,tVec] = OptimalTrajectoryDirectMat(N,timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,Ucart,mu,uMax)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method. 
% Method: direct transcription with variables in matricial form

% INPUT
% N: number of grid points
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% Ucart: control acceleration
% mu: gravitational parameter
% uMax: maximum control
%
% OUTPUT
% r: matrix with position
% v: matrix with velocity
% u: matrix with control acceleration
% tVec: vector with time 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters and variables for direct transcription

% Choose if to use a low order approximation for the collocation
% constraint or if to use the higher order Gauss-Lobatto scheme

LOW_ORDER = 1;

if LOW_ORDER
    
    % Find time step from number of nodes and time vector
    tVec = linspace(t0,tf,N);
    h = tVec(2)-tVec(1);
    
    % Write initial and final conditions
    inCond.r0 = r0;
    inCond.v0 = v0;
    inCond.t0 = t0;
    fCond.rf = rf;
    fCond.vf = vf;
    fCond.tf = tf;
    
    % Write initial guess vector: it must be similar in form to y that is found
    % by fmincon, so comprising the state and the control at each moment of time
    % y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
    % a 3x1 vector (velocity) and u is a 3x1 vector (control)
    % It is a matrix [3x3n x 1]
    
    % An interpolation is needed in order to have the values
    for i=1:3
        rInt(i,:) = spline(timeVec,R(i,:),tVec);
        vInt(i,:) = spline(timeVec,Vcart(i,:),tVec);
        uInt(i,:) = spline(timeVec,Ucart(i,:),tVec);
    end
    
    %Add some 'noise' to the initial guess data (a randomatic variation given
    %as a percentual of the original value)
%     a = -0.1;
%     b = 0.1;
%      for i=1:3
%       rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,N)).*rInt(i,:);
%       vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,N)).*vInt(i,:);
%       uInt(i,:) = uInt(i,:) + (a + (b-a).*rand(1,N)).*uInt(i,:);
%      end
    
    x0(1:3,1:N) = rInt;
    x0(4:6,1:N) = vInt;
    x0(7:9,1:N) = uInt;
    
    % Find the correct initial conditions for the lagrange multipliers from
    % fmincon
    options = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','TolX',1e-9,'TolCon',1e-8,'TolFun',1e-8,'MaxFunEvals',200000,'MaxIter',1000);
    xSol = fmincon(@(y) DirTranscFunctionFMINCONMatObj(y,N,tVec),x0,[],[],[],[],[],[],@(y) DirTranscFunctionFMINCONMatConst(y,inCond,fCond,N,h,mu,uMax),options);
    
    % Extract the obtained state and save it in matrices
    r = xSol(1:3,1:N);
    v = xSol(4:6,1:N);
    u = xSol(7:9,1:N);
    
else
    
    % Find time step from number of nodes and time vector
    tVec = linspace(t0,tf,2*N-1);
    % Time step (is between the two border nodes, not the middle node)
    h = tVec(3)-tVec(1);
    
    % Write initial and final conditions
    inCond.r0 = r0;
    inCond.v0 = v0;
    inCond.t0 = t0;
    fCond.rf = rf;
    fCond.vf = vf;
    fCond.tf = tf;
    
    % Write initial guess vector: it must be similar in form to y that is found
    % by fmincon, so comprising the state and the control at each moment of time
    % y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
    % a 3x1 vector (velocity) and u is a 3x1 vector (control)
    % It is a matrix [3x3n x 1]
    
    % An interpolation is needed in order to have the values
    for i=1:3
        rInt(i,:) = spline(timeVec,R(i,:),tVec);
        vInt(i,:) = spline(timeVec,Vcart(i,:),tVec);
        uInt(i,:) = spline(timeVec,Ucart(i,:),tVec);
    end
    
    %Add some 'noise' to the initial guess data (a randomatic variation given
    %as a percentual of the original value)
%     a = -0.5;
%     b = 0.5;
%      for i=1:3
%       rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,2*N-1)).*rInt(i,:);
%       vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,2*N-1)).*vInt(i,:);
%       uInt(i,:) = uInt(i,:) + (a + (b-a).*rand(1,2*N-1)).*uInt(i,:);
%      end
    
    x0(1:3,1:2*N-1) = rInt;
    x0(4:6,1:2*N-1) = vInt;
    x0(7:9,1:2*N-1) = uInt;
    
    % Find the correct initial conditions for the lagrange multipliers from
    % fmincon
    options = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','TolX',1e-8,'TolCon',1e-8,'TolFun',1e-8,'MaxFunEvals',200000,'MaxIter',100000);
    xSol = fmincon(@(y) DirTranscFunctionFMINCONGLObj(y,N,tVec),x0,[],[],[],[],[],[],@(y) DirTranscFunctionFMINCONGLConst(y,inCond,fCond,N,h,mu,uMax),options);
    
    % Extract the obtained state and save it in matrices
    r = xSol(1:3,1:2*N-1);
    v = xSol(4:6,1:2*N-1);
    u = xSol(7:9,1:2*N-1);
    
end

end