function [r,v,lambda,uVec,tVec] = OptimalTrajectoryIndirectFSOLVEMultipleLinUScaled(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,M,nT,mu,uMax)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method

% Method: indirect multiple shooting with FSOLVE, version for linear
% control formulation with scaled equations

% INPUT
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% lambdAd: lagrange multipliers
% M: number of grid points
% nT: number of integration points in a phase
% mu: gravitational parameter
%
% OUTPUT
% r: matrix with position
% v: matrix with velocity
% lambda: matrix with lagrange multipliers
% tVec: vector with time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define errors and tolerances for the shooting method

% Write initial and final conditions
inCond.r0 = r0;
inCond.v0 = v0;
inCond.t0 = t0;
fCond.rf = rf;
fCond.vf = vf;
fCond.tf = tf;
lambdaS = [norm(lambdaAd(1:3,1)) norm(lambdaAd(4:6,1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write initial guess of solution vector, which is made by M*4 vector
% variables (M*4*3 scalar variables)

% First interpolate to obtain the values at such time instants (with
% initial and final instants)
timeInt = linspace(t0,tf,M+1);

% An interpolation is needed in order to have the values
 for i=1:3
  rInt(i,:) = spline(timeVec,R(i,:),timeInt);
  vInt(i,:) = spline(timeVec,Vcart(i,:),timeInt);
  lambdarInt(i,:) = spline(timeVec,lambdaAd(i,:),timeInt);
  lambdavInt(i,:) = spline(timeVec,lambdaAd(i+3,:),timeInt);
 end

% %  Add some 'noise' to the initial guess data (a randomatic variation given
% % as a percentual of the original value)
% a = -0.001;
% b = 0.001;
%  for i=1:3
%   rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,M+1)).*rInt(i,:);
%   vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,M+1)).*vInt(i,:);
%   lambdarInt(i,:) = lambdarInt(i,:) + (a + (b-a).*rand(1,M+1)).*lambdarInt(i,:);
%   lambdavInt(i,:) = lambdavInt(i,:) + (a + (b-a).*rand(1,M+1)).*lambdavInt(i,:);
%  end
 
 % Write initial guess of solution vector (can choose if to pass or not
 % last element as variable -> select if 1:M or 1:M+1. Can also choose if
 % use the first element as opt variable or not -> select 2: instead of 1:)
 x0(1:3,:) = rInt(1:3,1:M);
 x0(4:6,:) = vInt(1:3,1:M);
 x0(7:9,:) = lambdarInt(1:3,1:M);
 x0(10:12,:) = lambdavInt(1:3,1:M);
 
% Pass the initial guess to fsolve
optionsFS = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter','TolFun',1e-3,'TolX',1e-2,'MaxFunEvals',5000,'MaxIter',4);
[xSol,fval] = fsolve(@(y) ShootingFunctionFSOLVEMultipleLinUScaled(y,inCond,fCond,mu,timeInt,M,nT,uMax,lambdaS),x0,optionsFS);

% Pass the initial guess to fmincon
%optionsFM = optimoptions('fmincon','Algorithm','interior-point','Display','iter','TolX',1e-8,'TolCon',1e-8,'TolFun',1e-8,'MaxFunEvals',100000,'MaxIter',100000);
%xSol = fmincon(@(y) ShootingFunctionFMINCONObj(y,M,timeInt,nT,mu),x0,[],[],[],[],[],[],@(y) ShootingFunctionFMINCONConst(y,inCond,fCond,mu,timeInt,M,nT),optionsFM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Obtain the solution (as matrix with initial values of r,v,lambda in each
% sub-interval)
rSol = xSol(1:3,:);
vSol = xSol(4:6,:);
lambdaSol(1:3,:) = xSol(7:9,:);
lambdaSol(4:6,:) = xSol(10:12,:);

% Version modified without r(1:3,1) as optimization variable (modify also
% the starting i in the cycle)
% i=1;
% [tOde,Y] = ode45(@(t,x) diffEq(t,x,mu),linspace(timeInt(i),timeInt(i+1),nT+1),[r0 v0 lambda0]);
% tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(1:end-1);
% r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,1:3);
% v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,4:6);
% lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(1:end-1,7:12);
      
% Execute M integrations on the different subintervals
for i=1:M 
    % Initial conditions for each sub-interval
    r_in = rSol(:,i);
    v_in = vSol(:,i);
    lambda_in = lambdaSol(:,i);
    
    % Integration on each sub-interval
    %[tOde,Y] = ode45(@(t,x) diffEq(t,x,mu),linspace(timeInt(i),timeInt(i+1),nT+1),[r_in;v_in;lambda_in]);
    timeInteg = linspace(timeInt(i),timeInt(i+1),nT+1);
    [tOde,Y] = rk4(@diffEqLinUScaled,timeInteg,xSol(:,i),1,mu,uMax,norm(r0),norm(v0),lambdaS(1),lambdaS(2),tf);
    Y = Y';
     
    % Saving of elements for each sub-interval (exclude the last point
    % because will be the first point of the next sub-interval integration) 
    tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(1:end-1);
    r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,1:3);
    v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,4:6);
    lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(1:end-1,7:12);
   
%     % Different way: saving of elements for each sub-interval (exclude the first point
%     % because will be the last point of the previous sub-interval integration) 
%     tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(2:end);
%     r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(2:end,1:3);
%     v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(2:end,4:6);
%     lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(2:end,7:12);
   
end

% Add the last element of the integration in the last sub-interval (because
% is not saved in the cycle), which corresponds to the last time point
tVec(end+1) = tOde(end);
r(end+1,1:3) = Y(end,1:3);
v(end+1,1:3) = Y(end,4:6);
lambda(end+1,1:6) = Y(end,7:12);

% %Add the first element of the integration in the first sub-interval (because
% %is not saved in the cycle), which corresponds to the first time point
% tVec = [t0 tVec];
% r = [rSol(1:3,1)';r];
% v = [vSol(1:3,1)';v];
% lambda = [lambdaSol(1:6,1)';lambda];

%--------------------------------------------------------------------------
% Find the control vector from lambda
%--------------------------------------------------------------------------
% Write lambdav
lambdav = lambda(:,4:6);

% Find the direction of the control from the adjoints
for i=1:length(tVec)
    uDir(i,:) = - lambdav(i,:)/norm(lambdav(i,:));
    
    % Calculate the Switching Function
    S = 1 - norm(lambdav(i,:));
    
    % Choose control modulus depending upon the Switching Function
    if S < 0
        uMod = uMax;
    else
        uMod = 0;
    end
    
    % Define the vector of control
    uVec(i,:) = uMod*uDir(i,:);
    
end

end

