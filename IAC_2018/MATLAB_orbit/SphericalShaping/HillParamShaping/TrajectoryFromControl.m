function [R_ver,V_ver] = TrajectoryFromControl(tVec,r0,v0,Ucart,mu,r,v)
%
% Function TrajectoryFromControl: integrate the dynamics equation of a low
% thrust spacecraft with a known control profile in time
% The integration is done in cartesian coordinates, as well as the inputs
% and outputs
%
% INPUT
% tVec: time vector
% r0,v0: initial values for position and velocity
% Ucart: control profile matrix (cartesian coordinates) of size (3 x n)
% mu: gravitational parameter
% r,v: matrices with the original position and velocities (size (3 x n))
%
% OUTPUT
% R_ver,V_ver: matrix of cartesian position and velocity during time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate the trajectory with the control

% Transform into column vector initial conditions for requirement of ode45
r0 = r0(:);
v0 = v0(:);

% Integrate with ode45
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'NormControl','on','MaxStep',1);
%[tVec,Y] = ode45(@(t,x) diffEqWithControl(t,x,mu,Ucart,tVec),tVec,[r0;v0],options);

% %Integrate with rk4
% [tVec,Y] = rk4(@diffEqWithControl,tVec,[r0;v0],1,mu,Ucart',tVec);
% Y = Y';

% Extract position and velocity
% R_ver = (Y(:,1:3))';
% V_ver = (Y(:,4:6))';

% Integrate to obtain sol
sol = ode113(@(t,x) diffEqWithControl(t,x,mu,Ucart,tVec),tVec,[r0;v0],options);
[Y,YP] = deval(sol,tVec);
R_ver = (Y(1:3,:));
V_ver = (Y(4:6,:));

% Plot results obtained and compare with input vectors (r,v) that can be
% from an inverse method or from an optimal control method
%Radius vector
figure;
subplot(3,1,1);
plot(tVec,r(1,:),'g','linewidth',2);
hold on;
plot(tVec,R_ver(1,:),'r','linewidth',2);
ylabel('R x');
legend('Pseudo Equin Shaping','Int From Control','location','northeast');
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
legend('Pseudo Equin Shaping','Int From Control','location','northeast');
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

%Control acceleration recovered from integration
Acc_ver = YP(4:6,:);
RMod_ver = sqrt(R_ver(1,:).^2  + R_ver(2,:).^2 + R_ver(3,:).^2);
for i=1:3
   ACont(i,:) = Acc_ver(i,:) + (mu./RMod_ver.^3).*R_ver(i,:);
end

% Phisical acceleration recovered from original matrices of state and
% control
RMod = sqrt(r(1,:).^2  + r(2,:).^2 + r(3,:).^2);
for i=1:3
   Acc(i,:) = Ucart(i,:) - (mu./RMod.^3).*r(i,:);
end

% Plot control acceleration of original trajectory and of verified
% trajectory
figure;
subplot(3,1,1);
plot(tVec,Ucart(1,:),'g','linewidth',2);
hold on;
plot(tVec,ACont(1,:),'r','linewidth',2);
ylabel('U x');
legend('Pseudo Equin Shaping','Int From Control','location','northeast');
subplot(3,1,2);
plot(tVec,Ucart(2,:),'g','linewidth',2);
hold on;
plot(tVec,ACont(2,:),'r','linewidth',2);
ylabel('U y');
subplot(3,1,3);
plot(tVec,Ucart(3,:),'g','linewidth',2);
hold on;
plot(tVec,ACont(3,:),'r','linewidth',2);
ylabel('U z');
suptitle('Verification for Control Acceleration vector');
xlabel('Time');

% Plot phisical acceleration of original trajectory and of verified
% trajectory
figure;
subplot(3,1,1);
plot(tVec,Acc(1,:),'g','linewidth',2);
hold on;
plot(tVec,Acc_ver(1,:),'r','linewidth',2);
ylabel('U x');
legend('Pseudo Equin Shaping','Int From Control','location','northeast');
subplot(3,1,2);
plot(tVec,Acc(2,:),'g','linewidth',2);
hold on;
plot(tVec,Acc_ver(2,:),'r','linewidth',2);
ylabel('U y');
subplot(3,1,3);
plot(tVec,Acc(3,:),'g','linewidth',2);
hold on;
plot(tVec,Acc_ver(3,:),'r','linewidth',2);
ylabel('U z');
suptitle('Verification for Acceleration vector');
xlabel('Time');
end

