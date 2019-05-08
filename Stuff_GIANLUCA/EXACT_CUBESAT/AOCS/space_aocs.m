function [M,P,info] = space_aocs(x,ep)

%% space_aocs: Attitude and orbit control system model
%
% [M,P,varargout] = space_aocs(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = pointing accuracy, deg
%       * x(2) = offset between centre of gravity and centre of pressure, m
%       * x(3) = area normal to velocity vector, m^2
%       * x(4) = area normal to sun line, m^2
%       * x(5) = reflectance factor
%       * x(6) = spacecraft residual dipole, A*m^2
%       * x(7) = drag coefficient
%       * x(8) = actuator type: [0 0.5) = magtorq, [0.5 1] = thruster
%       * x(9) = specific impulse, s
%       * x(10) = burn time for momentum dumping, s
%       * x(11) = slew angle, deg
%       * x(12) = time for slew maneuvre, s
%       * x(13) = initial tumbling rate, rad/s
%       * x(14) = spin angle, deg      
% * ep: Environmental parameters
%       * ep(1) = orbit period, s
%       * ep(2) = mean planet's magnetic field stregth, T
%       * ep(3) = mean incident solar radiation, W/m^2
%       * ep(4) = mean dynamic pressure, Pa
%       * ep(5) = mean gravitational field strength, deg
%       * ep(6) = angle between spacecraft z axis and nadir vector, deg
%       * ep(7) = spacecraft moment of inertia about x axis, kg*m^2
%       * ep(8) = spacecraft moment of inertia about y axis, kg*m^2
%       * ep(9) = spacecraft moment of inertia about z axis, kg*m^2
%       * ep(10) = thruster moment arm, m
%       * ep(11) = number of thrusters
%       * ep(12) = number of dumping maneuvres during S/C lifetime
%       * ep(13) = number of slew maneuvres during S/C lifetime
%
%%  Output:
% * M = total mass of the aocs subsystem
% * P = power consumption of the aocs subsystem
% * info:
%       * info.M_propellant = propellant mass 
%       * info.max_F_rct = maximum force per thruster in RCS
%
%% Author: Simone Alicino, 2013

% Design variables
accuracy = x(1);
lsp = x(2);
Av = x(3);
As = x(4);
q = x(5);
m = x(6);
Cd = x(7);
if x(8) < 0.5
    actuator = 'magtorquer';
else
    actuator = 'thruster';
end
Isp = x(9);
tdump = x(10);
slew = x(11)*pi/180;
t_slew = x(12);
spin_rate = x(13);
spin_angle = x(14)*pi/180;

% Environmental parameters
period = ep(1);
B = ep(2);
Is = ep(3);
dyn_pres = ep(4);
G = ep(5);
theta = ep(6);
Ix = ep(7);
Iy = ep(8);
Iz = ep(9);
L = ep(10);
n_thrusters = ep(11);
ndump = ep(12);
nslew = ep(13);

% Hardcoded parameters
c = 299792458;  % speed of light, m/s
M_propellant = 0;
F_rct = [];

% Solar pressure
Ts = Is*As*lsp*(1 + q)/c;

% Magnetic torque
Tm = m*B;

% Gravity-gradient torque
Tg = G*abs(Iz - min([Ix Iy]))*sind(2*theta)/2;

% Aerodynamic drag torque
Ta = dyn_pres*Cd*Av*lsp;

% Total disturbance torque
Td = Ts + Tm + Tg + Ta;

% Momentum storage
Hd = Td*period/(4*accuracy*pi/180);
[M_wheel(1), P_wheel(1)] = wheels(Hd);

% Momentum dumping
switch actuator
    case 'magtorquer'
        D_dumping = Td/B;
        [M_magtorq(1), P_magtorq(1)] = magnetorquers(D_dumping);
    case 'thruster'
        F_dumping = Hd/(tdump*L);
        [M_thruster(1), P_thruster(1)] = thrusters(F_dumping,tdump,Isp,ndump);
end
if strcmpi(actuator,'thruster')
    M_propellant = M_propellant + M_thruster(1);
    F_rct = [F_rct, F_dumping];
end

% De-tumbling maneuver
H_tumbl = Iz*spin_rate;
T_tumbl = H_tumbl*spin_rate/(2*spin_angle);
[M_tumbl(1), P_tumbl(1)] = wheels(H_tumbl);
switch actuator
    case 'magtorquer'
        D_tumbl = T_tumbl/B/cosd(45);
        [M_tumbl(2), P_tumbl(2)] = magnetorquers(D_tumbl);
    case 'thruster'
%         F_tumbl = T_tumbl/(n_thrusters*L);
        F_tumbl = T_tumbl/L;
        [M_tumbl(2), P_tumbl(2)] = thrusters(F_tumbl,0.5*t_slew,Isp,nslew);
end
[Mtumbl,i] = min(M_tumbl);
Ptumbl = P_tumbl(i);
if strcmpi(actuator,'thruster') && i == 2
    M_propellant = M_propellant + M_tumbl(2);
    F_rct = [F_rct, F_tumbl];
end

% Slew maneuver
slew_rate = 2*slew/t_slew;
T_slew = 4*slew*Iz/t_slew^2;
H_slew = T_slew*t_slew;
[M_slew(1), P_slew(1)] = wheels(H_slew);
switch actuator
    case 'magtorquer'
        D_slew = T_slew/B;
        [M_slew(2), P_slew(2)] = magnetorquers(D_slew);
    case 'thruster'
        F_slew = Iz*slew_rate/L;
        [M_slew(2), P_slew(2)] = thrusters(F_slew,0.5*t_slew,Isp,nslew);
end
[Mslew,i] = min(M_slew);
Pslew = P_slew(i);
if strcmpi(actuator,'thruster') && i == 2
    M_propellant = M_propellant + M_slew(2);
    F_rct = [F_rct, F_slew];
end

Mact = [Mtumbl,Mslew];
Pact = [Ptumbl,Pslew];
[M_actuator,i] = max(Mact);
P_actuator = Pact(i);

% switch actuator
%     case 'magtorquer'
%         M_actuator = max(M_magtorq);
%         P_actuator = max(P_magtorq);
%     case 'thruster'
%         M_actuator = sum(M_thruster);
%         P_actuator = sum(P_thruster);
% end

% Total mass and power
M = M_wheel + M_actuator;
P = P_wheel + P_actuator;
info.M_propellant = M_propellant;
F_rct(isempty(F_rct)) = 0;
max_F_rct = max(F_rct)/n_thrusters;
max_F_rct(isnan(max_F_rct)) = 0;
info.max_F_rct = max_F_rct;

end

function [M_wheel, P_wheel, wheel] = wheels(H)

M_wheel = interp1([0.0016 0.4 400],[0.072 2 20],H,'linear','extrap');
P_wheel = interp1([0.0016 0.4 400],[0.465 10 110],H,'linear','extrap');

wheel.mass = M_wheel;
wheel.power = P_wheel;
wheel.momentum = H;

end

function [M_rct, P_rct, rct] = thrusters(F,t,Isp,n)

% Propellant mass only!
% F = thrust, N
% t = maneuvre time, s
% Isp = specific impulse, s
% n = number of maneuvres in spacecraft lifetime

g0 = 9.80665;
Mp = F*t/(Isp*g0)*n;

M_rct = Mp;
P_rct = 0;

rct.mass = M_rct;
rct.power = P_rct;

end

function [M_magtorq, P_magtorq, magtorq] = magnetorquers(D)

M_magtorq = interp1([0.06 1 4000],[0.0835 0.4 50],D,'linear','extrap');
P_magtorq = interp1([0.06 1 4000],[0.155 0.6 16],D,'linear','extrap');

magtorq.mass = M_magtorq;
magtorq.power = P_magtorq;
magtorq.dipole = D;

end
