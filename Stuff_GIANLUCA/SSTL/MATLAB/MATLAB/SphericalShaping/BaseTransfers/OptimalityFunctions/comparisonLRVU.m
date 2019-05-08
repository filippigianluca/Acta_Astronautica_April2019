function comparisonLRVU(timeINV,timeOPT,rINV,rOPT,vINV,vOPT,uINV,uOPT,varargin)
%
%comparisonLRVU: function that plot the adjoints,radius,velocity and
%control of inverse method and of optimal control
%
% INPUT
% vectors and matrices corresponding to inverse method and optimal control
% varargin because for indirect method also the adjoints are plotted
%
%
% OUTPUT:
% no output, just graph plotting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

% Set constants
AU2km = 149597870.66;
T_sid = 86164.1004;

% Radius vector plot
figure;
subplot(3,1,1);
plot(timeINV,rINV(1,:),'b','linewidth',2);
hold on;
plot(timeOPT,rOPT(1,:),'g','linewidth',2);
ylabel('R x');
legend('Spherical Shaping','Optimal Control','location','northeast');
subplot(3,1,2);
plot(timeINV,rINV(2,:),'b','linewidth',2);
hold on;
plot(timeOPT,rOPT(2,:),'g','linewidth',2);
ylabel('R y');
subplot(3,1,3);
plot(timeINV,rINV(3,:),'b','linewidth',2);
hold on;
plot(timeOPT,rOPT(3,:),'g','linewidth',2);
ylabel('R z');
suptitle('Radius vector comparison Spherical Shaping-Optimal Control');
xlabel('Time');

% Velocity vector plot
figure;
subplot(3,1,1);
plot(timeINV,vINV(1,:),'b','linewidth',2);
hold on;
plot(timeOPT,vOPT(1,:),'g','linewidth',2);
ylabel('V x');
legend('Spherical Shaping','Optimal Control','location','northeast');
subplot(3,1,2);
plot(timeINV,vINV(2,:),'b','linewidth',2);
hold on;
plot(timeOPT,vOPT(2,:),'g','linewidth',2);
ylabel('V y');
subplot(3,1,3);
plot(timeINV,vINV(3,:),'b','linewidth',2);
hold on;
plot(timeOPT,vOPT(3,:),'g','linewidth',2);
ylabel('V z');
suptitle('Velocity vector comparison Spherical Shaping-Optimal Control');
xlabel('Time');

% Control vector plot
figure;
subplot(3,1,1);
plot(timeINV,uINV(1,:)*AU2km/T_sid^2,'b','linewidth',2);
hold on;
plot(timeOPT,uOPT(1,:)*AU2km/T_sid^2,'g','linewidth',2);
ylabel('U x  [km/s^2]');
legend('Spherical Shaping','Optimal Control','location','northeast');
subplot(3,1,2);
plot(timeINV,uINV(2,:)*AU2km/T_sid^2,'b','linewidth',2);
hold on;
plot(timeOPT,uOPT(2,:)*AU2km/T_sid^2,'g','linewidth',2);
ylabel('U y  [km/s^2]');
subplot(3,1,3);
plot(timeINV,uINV(3,:)*AU2km/T_sid^2,'b','linewidth',2);
hold on;
plot(timeOPT,uOPT(3,:)*AU2km/T_sid^2,'g','linewidth',2);
ylabel('U z  [km/s^2]');
suptitle('Control vector comparison Spherical Shaping-Optimal Control');
xlabel('Time');

% Plot the modulus of the control acceleration
figure;
plot(timeINV,sqrt(uINV(1,:).^2 + uINV(2,:).^2 + uINV(3,:).^2)*AU2km/T_sid^2,'b','LineWidth',2.0);
hold on;
plot(timeOPT,sqrt(uOPT(1,:).^2 + uOPT(2,:).^2 + uOPT(3,:).^2)*AU2km/T_sid^2,'g','LineWidth',2.0);
ylabel('U [km/s^2]');
xlabel('Time [days]');
legend('Spherical Shaping','Optimal Control');
suptitle('Control modulus comparison Spherical Shaping-Optimal Control');

% Check number of varaibles in input
if (nargin == 10)
    
    % Assign variables from varargin
    lINV = varargin{1};
    lOPT = varargin{2};
    
    % Adjoints vector plot
    figure;
    subplot(3,1,1);
    plot(timeINV,lINV(1,:),'b','linewidth',2);
    hold on;
    plot(timeOPT,lOPT(1,:),'g','linewidth',2);
    ylabel('\lambda rx');
    legend('Spherical Shaping','Optimal Control','location','northeast');
    subplot(3,1,2);
    plot(timeINV,lINV(2,:),'b','linewidth',2);
    hold on;
    plot(timeOPT,lOPT(2,:),'g','linewidth',2);
    ylabel('\lambda ry');
    subplot(3,1,3);
    plot(timeINV,lINV(3,:),'b','linewidth',2);
    hold on;
    plot(timeOPT,lOPT(3,:),'g','linewidth',2);
    ylabel('\lambda rz');
    suptitle('r adjoints comparison Spherical Shaping-Optimal Control');
    xlabel('Time');
end

end

