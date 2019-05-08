function [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Ucart,mu,r,v)
%
% Function TrajectoryFromControl: integrate the dynamics equation of a low
% thrust spacecraft with a known control profile in time
% The integration is done in cartesian coordinates, as well as the inputs
% and outputs
%
% INPUT
% time: time vector
% r0,v0: initial values for position and velocity
% Ucart: control profile matrix (cartesian coordinates) of size (3 x n)
% mu: gravitational parameter
%
% OUTPUT
% R,V: matrix of cartesian position and velocity during time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate the trajectory with the control

% Transform into column vector initial conditions for requirement of ode45
r0 = r0(:);
v0 = v0(:);

% Integrate with ode45
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[tVec,Y] = ode45(@(t,x) diffEqWithControl(t,x,mu,Ucart,tVec),tVec,[r0;v0],options);

% %Integrate with rk4
%[tVec,Y] = rk4(@diffEqWithControl,tVec,[r0;v0],1,mu,Ucart,tVec);
%Y = Y';

% Extract position and velocity
R_ver = (Y(:,1:3))';
V_ver = (Y(:,4:6))';

% Plot results obtained and compare with input vectors (r,v) that can be
% from an inverse method or from an optimal control method
%Radius vector
figure;
subplot(3,1,1);
plot(tVec,r(1,:),'g','linewidth',2);
hold on;
plot(tVec,R_ver(1,:),'r','linewidth',2);
ylabel('R x');
legend('Optimal Control Result','Int From Control','location','northeast');
subplot(3,1,2);
plot(tVec,r(2,:),'g','linewidth',2);
hold on;
plot(tVec,R_ver(2,:),'r','linewidth',2);
ylabel('R y');
subplot(3,1,3);
plot(tVec,r(3,:),'g','linewidth',2);
hold on;
plot(tVec,R_ver(3,:),'r','linewidth',2);
ylabel('R z');
suptitle('Verification for Radius vector');
xlabel('Time');

%Velocity vector
figure;
subplot(3,1,1);
plot(tVec,v(1,:),'g','linewidth',2);
hold on;
plot(tVec,V_ver(1,:),'r','linewidth',2);
ylabel('V x');
legend('Optimal Control Result','Int From Control','location','northeast');
subplot(3,1,2);
plot(tVec,v(2,:),'g','linewidth',2);
hold on;
plot(tVec,V_ver(2,:),'r','linewidth',2);
ylabel('V y');
subplot(3,1,3);
plot(tVec,v(3,:),'g','linewidth',2);
hold on;
plot(tVec,V_ver(3,:),'r','linewidth',2);
ylabel('V z');
suptitle('Verification for Velocity vector');
xlabel('Time');

end

