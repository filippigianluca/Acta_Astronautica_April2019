function [M,P,info] = space_thermal(x,ep)

%% space_thermal: Thermal model
%
%   [M,P,info] = space_thermal(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = absorptivity of radiator
%       * x(2) = emissivity of radiator
%       * x(3) = absorptivity of MLI
%       * x(4) = emissivity of MLI
%       * x(5) = efficiency of MLI
%       * x(6) = absorptivity degradation (to compute EOL alpha), %/yr
%       * x(7) = specific mass of radiator, kg/m^2
%       * x(8) = specific mass of MLI, kg/m^2
%       * x(9) = emissivity of louvers
%       * x(10) = specific mass of louvers, kg/m^2
%       * x(11) = specific mass of paint, kg/m^2
%       * x(12) = specific mass of heater, kg/W
%       * x(13) = length of thermal link, m
%       * x(14) = linear density of thermal link device, kg/m
%
% * ep: Environmental parameters
%       * ep(1:2) = minimum:maximum internal heat dissipation, W
%       * ep(3:4) = minimum:maximum temperature limit, C
%       * ep(5:6) = minimun:maximum IR emission of the planet, W/m^2
%       * ep(7:8) = minimum:maximum albedo
%       * ep(9:10) = minimum:maximum Solar intensity, W/m^2
%       * ep(11:12) = minimum:maximum planet view factor
%       * ep(13:14) = minimum:maximum albedo view factor
%       * ep(15:16) = minimum:maximum sun view factor
%       * ep(17) = total area of spacecraft, m^2
%       * ep(18) = mission duration, years
%
%% Outputs:
% * M = total mass of the thermal subsystem
% * P = electrical power required by the thermal subsystem
% * info:
%       * info.mass = mass inventory
%       * info.temperature = temperatures
%       * info.area = area
%
%% Author: Simone Alicino, 2013

% Design parameters
alpha_r = x(1);
epsilon_r = x(2);
alpha_m = x(3);
epsilon_m = x(4);
eta_mli = x(5);
degradation = x(6);
rho_radiator = x(7);
rho_mli = x(8);
epsilon_l = x(9);
rho_louvers = x(10);
rho_paint = x(11);
k_heater = x(12);
l = x(13);
rho_link = x(14);

% Environmental parameters
Q_cold = ep(1);
Q_hot = ep(2);
T_min = 273.15 + ep(3);
T_max = 273.15 + ep(4);
IR_planet = ep(5:6);
a = ep(7:8);
G = ep(9:10);
F_planet = ep(11:12);
F_albedo = ep(13:14);
F_sun = ep(15:16);
A_spacecraft = ep(17);
life = ep(18);

% Hardcoded parameters
sigma = 5.67e-8;

% Solar, planet, albedo fluxes
IS = G;
IP = IR_planet;
IA = a.*IS;
Q_sun = IS.*F_sun;
Q_planet = IP.*F_planet;
Q_albedo = IA.*F_albedo;

% Worst case hot and cold heat fluxes
Q_sun_max = max(Q_sun);
Q_albedo_max = max(Q_albedo);
Q_planet_max = max(Q_planet);

Q_sun_min = min(Q_sun);
Q_albedo_min = min(Q_albedo);
Q_planet_min = min(Q_planet);

% Effective MLI absorptivity and emissivity
alpha_eff = alpha_m*eta_mli/(eta_mli + epsilon_m);
epsilon_eff = epsilon_m*eta_mli/(eta_mli + epsilon_m);

% HOT CASE: Radiator sizing

% EOL absorptivity
alpha_r_EOL = alpha_r*(1 + degradation*life/100);
alpha_m_EOL = alpha_m;
% alpha_m_EOL = alpha_m*(1 + degradation*life/100);
alpha_eff_EOL = alpha_m_EOL*eta_mli/(eta_mli + epsilon_m);

Dalpha = alpha_r_EOL - alpha_eff_EOL;
Depsilon = epsilon_r - epsilon_eff;

% Radiator
Q_alpha = Q_sun_max + Q_albedo_max;
Q_epsilon = Q_planet_max - sigma*T_max^4;
A_radiator = ( A_spacecraft*alpha_eff_EOL*Q_alpha + A_spacecraft*epsilon_eff*Q_epsilon + Q_hot)/...
    ( -Depsilon*Q_epsilon - Dalpha*Q_alpha );
A_radiator(A_radiator < 0) = 0;
A_deployed_radiator = 0;
if A_radiator > A_spacecraft
    A_deployed_radiator = A_radiator - A_spacecraft;
    A_radiator = A_spacecraft;
end
M_body_radiator = A_radiator*rho_radiator;
M_deployed_radiator = A_deployed_radiator*rho_radiator;
M_radiator = M_body_radiator + M_deployed_radiator;

% COLD CASE: heater sizing

% Equivalent absorptivity and emissivity
alpha1 = alpha_eff + A_radiator/A_spacecraft*(alpha_r - alpha_eff);
epsilon1 = epsilon_eff + A_radiator/A_spacecraft*(epsilon_r - epsilon_eff);

% Louvers
M_louvers = 0;
X_louvers = (alpha1*(Q_sun_min + Q_albedo_min) + epsilon1*Q_planet_min + Q_cold/A_spacecraft)/(sigma*epsilon1);
T_radiator_cold = nthroot(X_louvers, 4);
if T_radiator_cold < T_min
    M_louvers = A_radiator*rho_louvers;
    epsilon_r = epsilon_l;
    epsilon1 = epsilon_eff + A_radiator/A_spacecraft*(epsilon_r - epsilon_eff);
%     alpha_r = alpha_r;
%     alpha1 = alpha_eff + A_radiatior/A_spacecraft*(alpha_r - alpha_eff);
end

% Heater
P_heater = A_spacecraft*( epsilon1*(sigma*T_min^4 - Q_planet_min) - alpha1*(Q_sun_min + Q_albedo_min) ) - Q_cold;
P_heater(P_heater < 0) = 0;
M_heater = P_heater*k_heater;

% OTHER HARDWARE

% Multilayer Insulation (MLI)
A_mli = A_spacecraft - A_radiator;
M_mli = A_mli*rho_mli;

% Paint (interior)
M_paint = A_spacecraft*rho_paint;

% Thermal link
M_link = l*rho_link;

% Total mass and power consumption
M = M_radiator + M_heater + M_mli + M_louvers + M_link + M_paint;
P = P_heater;

mass.radiator = M_body_radiator;
mass.deployed_radiator = M_deployed_radiator;
mass.louvers = M_louvers;
mass.heater = M_heater;
mass.mli = M_mli;
mass.paint = M_paint;
mass.link = M_link;
mass.total = M;
info.mass= mass;

area.radiator = A_radiator;
area.deployed_radiator = A_deployed_radiator;
area.mli = A_mli;
info.area = area;

temperature.upper_limit = T_max - 273.15;
temperature.lower_limit = T_min - 273.15;
info.temperature = temperature;

end
