% Test_Shape_1: test on the shape for Hill parameters (plus verifiation of
% shape for pseudo-equinoctial)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and parameters
clc;
clear all;
close all;

% Load the data
LoadData;

% Set the parameters of the spacecraft
m_i = 1000;                  % kg
Ceff = 3000*9.81*1e-3;       % Ceff = Isp*g0 [km/s]

% Set tolerance,step and max iterations
toll = 1e-3;                % at least 1e-3 (1 day of error for TOF=100d)
thetaStep = 2*pi/200;
maxIter = 600;

% Set the departure time, TOF (MJD2000 and days), and num revs
t_dep = 8000;
TOF = 900;
timeStep = 1;
nr = 1;

% Estimate a maximum thrust limit
DV = 10*km2AU*T_sid;                 %[AU/day]
uMax = DV/TOF;                       %[AU/day^2]

%-------------------------------------------------------------------------
% Departure and Arrival bodies
%-------------------------------------------------------------------------
% Define Earth as departure body
DepBody = body(1);

% Select the target body
choice = menu('Select the target body','Mars','1989ML','Tempel-1','Neptune');
switch choice
    case 1
       ArrBody = body(2); 
    case 2
       ArrBody = body(3);
    case 3
       ArrBody = body(4); 
    case 4
       ArrBody = body(5); 
end


%% Integrate the trajectory with control

% Choose if constant thrust of thrust diminishing with inverse of distance
% squared from Sun
CONST_THRUST = 1;

% Initial conditions (from orbital elements)
% Find position and velocity of Departure Body
M_dep = DepBody.M0 + DepBody.n*(t_dep-0);
Theta_dep = M2theta(M_dep,DepBody.e);

% Transform the orbital elements into cartesian coordinates
[dep_rC,dep_vC] = KeplElem2rv(DepBody.a,DepBody.e,DepBody.i,DepBody.w,DepBody.Omega,Theta_dep,mu_S);
   
% Set initial conditions
a0 = DepBody.a;
x0 = [dep_rC;dep_vC];
r0 = norm(dep_rC);
nStep = TOF;
tVec = linspace(t_dep,t_dep + TOF,nStep);
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

%-------------------------------------------------------------------------
% Choose which integration to perform
%-------------------------------------------------------------------------  

if CONST_THRUST
    %-------------------------------------------------------------------------
    % Constant thrust
    %-------------------------------------------------------------------------
    [tOde,YOde] = ode45(@(t,x) diffEqTangU(t,x,mu_S,uMax), tVec, x0);
    rOde = YOde(:,1:3);
    vOde = YOde(:,4:6);
    
    % Extract the Keplerian elements
    for i=1:nStep
        Cart = [rOde(i,:) vOde(i,:)];
        [Kep] = Cart2Kep(Cart,mu_S);
        Pseudo(:,i) = Kep2Pseudo(Kep);
        Hill(:,i) = Kep2Hill(Kep,mu_S);
    end
    
    % Find the control
    for i=1:nStep
        normV = norm(vOde(i,:));
        unitV = vOde(i,:)/normV;
        Uopt(i,:) = uMax*unitV;
    end
    
    % Plot trajectory obtained
    plotTrajectoryCartesianState(DepBody,ArrBody,t_dep,TOF,rOde,Uopt');
    
else
    %-------------------------------------------------------------------------
    % Thrust decreasing with inverse of distance squared
    %-------------------------------------------------------------------------
    [tOde,YOde] = ode45(@(t,x) diffEqTangUDec(t,x,mu_S,uMax,r0), tVec, x0);
    rOde = YOde(:,1:3);
    vOde = YOde(:,4:6);
    
    % Extract the Keplerian elements
    for i=1:nStep
        Cart = [rOde(i,:) vOde(i,:)];
        [Kep] = Cart2Kep(Cart,mu_S);
        Pseudo(:,i) = Kep2Pseudo(Kep);
        Hill(:,i) = Kep2Hill(Kep,mu_S);
        %XParam(:,i) = Kep2XParam(Kep,mu_S,a0);
    end
    
    % Find the control
    for i=1:nStep
        modV = norm(vOde(i,:));
        unitV = vOde(i,:)/modV;
        modR = norm(rOde(i,:));
        Uopt(i,:) = (uMax*r0^2/modR^2)*unitV;
    end
    
    % Plot trajectory obtained
    plotTrajectoryCartesianState(DepBody,ArrBody,t_dep,TOF,rOde,Uopt');
    
end

%% Convert L such that is always increasing
% The algorithm takes into account that the numerical values of L are
% between [0-2pi]
L = Pseudo(6,:);

cont = 0;
toll = pi/5;
for i=1:nStep-1
   % Find discontinuities
   if ( ( (L(i) > 2*pi - toll) && (L(i) < 2*pi)) && ( (L(i+1) > 0) && (L(i+1) < toll)) )
    index(cont+1) = i+1;
    cont = cont + 1;
   end 
end

for i=1:length(index)
   % Add 2pi to the correspondent part
   L(index(i):end) = L(index(i):end) + 2*pi;
end

%% %% Convert u such that is always increasing
% The algorithm takes into account that the numerical values of u are
% between [0-2pi]
u = Hill(4,:);

cont = 0;
toll = pi/5;
for i=1:nStep-1
   % Find discontinuities
   if ( ( (u(i) > 2*pi - toll) && (u(i) < 2*pi)) && ( (u(i+1) > 0) && (u(i+1) < toll)) )
    index(cont+1) = i+1;
    cont = cont + 1;
   end 
end

for i=1:length(index)
   % Add 2pi to the correspondent part
   u(index(i):end) = u(index(i):end) + 2*pi;
end

%% Set the parameters vectors for the shape (with random values) for Pseudo equin
alpha0 = Pseudo(1:5,1);
alphaf = Pseudo(1:5,end);
lambda = [0.03 0.05 0.05 0.02 0.02];
L0 = L(1);

if CONST_THRUST
    for i=1:5
        alpha(i,:) = alpha0(i) + alphaf(i).*exp(lambda(i).*(L - L0));
    end
else
    for i=1:5
        alpha(i,:) = alpha0(i) + alphaf(i)*(L - L0) + lambda(i)*sin(L-L0);
    end
end

%% %% Set the parameters vectors for the shape (with random values) for Hill elem
beta0 = Hill(1:6,1);
betaf = Hill(1:6,end);
lambda = [0.03 0.05 0.05 0.02 0.02 0.01];
u0 = u(1);

if CONST_THRUST
    for i=1:6
        beta(i,:) = beta0(i) + betaf(i).*exp(lambda(i).*(u - u0));
    end
else
    for i=1:6
        beta(i,:) = beta0(i) + betaf(i)*(u - u0) + lambda(i)*sin(u-u0);
    end
end

%% Plot the results

%--------------------------------------------------------------------------
% Pseudo equinoctial elements plotting
%--------------------------------------------------------------------------
% First three
figure;
subplot(3,1,1);
plot(tOde,Pseudo(1,:),'b','linewidth',2);
hold on;
plot(tOde,alpha(1,:),'g','linewidth',2);
ylabel('p');
legend('Pseudo from Int','Shaped Pseudo','location','northeast');
subplot(3,1,2);
plot(tOde,Pseudo(2,:),'b','linewidth',2);
hold on;
plot(tOde,alpha(2,:),'g','linewidth',2);
ylabel('f');
subplot(3,1,3);
plot(tOde,Pseudo(3,:),'b','linewidth',2);
hold on;
plot(tOde,alpha(3,:),'g','linewidth',2);
ylabel('g');
suptitle('Comparison Pseudo elem from integration vs Pseudo elem from shaping');
xlabel('Time');
        
% Last three
figure;
subplot(3,1,1);
plot(tOde,Pseudo(4,:),'b','linewidth',2);
hold on;
plot(tOde,alpha(4,:),'g','linewidth',2);
ylabel('h');
legend('Pseudo from Int','Shaped Pseudo','location','northeast');
subplot(3,1,2);
plot(tOde,Pseudo(5,:),'b','linewidth',2);
hold on;
plot(tOde,alpha(5,:),'g','linewidth',2);
ylabel('k');
subplot(3,1,3);
plot(tOde,L,'b','linewidth',2);
ylabel('L');
suptitle('Comparison Pseudo elem from integration vs Pseudo elem from shaping');
xlabel('Time');
        
%--------------------------------------------------------------------------
% Hill elements plotting
%--------------------------------------------------------------------------
% First three
figure;
subplot(3,1,1);
plot(tOde,Hill(1,:),'b','linewidth',2);
hold on;
plot(tOde,beta(1,:),'g','linewidth',2);
ylabel('r');
legend('Hill from Int','Shaped Hill','location','northeast');
subplot(3,1,2);
plot(tOde,Hill(2,:),'b','linewidth',2);
hold on;
plot(tOde,beta(2,:),'g','linewidth',2);
ylabel('p');
subplot(3,1,3);
plot(tOde,Hill(3,:),'b','linewidth',2);
hold on;
plot(tOde,beta(3,:),'g','linewidth',2);
ylabel('G');
suptitle('Comparison Hill elem from integration vs Hill elem from shaping');
xlabel('Time');
        
% Last three
figure;
subplot(3,1,1);
plot(tOde,u,'b','linewidth',2);
hold on;
ylabel('u');
legend('Hill from Int','Shaped Hill','location','northeast');
subplot(3,1,2);
plot(tOde,Hill(5,:),'b','linewidth',2);
hold on;
plot(tOde,beta(5,:),'g','linewidth',2);
ylabel('H');
subplot(3,1,3);
plot(tOde,Hill(6,:),'b','linewidth',2);
hold on;
plot(tOde,beta(6,:),'g','linewidth',2);
ylabel('h');
suptitle('Comparison Hill elem from integration vs Hill elem from shaping');
xlabel('Time');
