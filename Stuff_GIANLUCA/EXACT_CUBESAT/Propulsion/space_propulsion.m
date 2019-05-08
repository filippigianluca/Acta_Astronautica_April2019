function [M,P,info] = space_propulsion(x,ep)

%% space_propulsion: Chemical propulsion model
%
%   [M,P,info] = space_propulsion(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = total initial spacecraft mass, kg
%       * x(2) = specific impluse, s [200 350]
%       * x(3) = combustion chamber pressure, Pa
%       * x(4) = blowdown ratio
%       * x(5) = exit/throat area ratio (epsilon)
%       * x(6) = tank material: [0 0.33] aluminum, (0.33 0.66) steel, [0.66 1] titanium
%       * x(7) = uncertainty on material properties, %
%       * x(8) = pressurant gas: [0 0.5) helium, [0.5 1] nitrogen
%       * x(9) = initial pressurant pressure, Pa
%       * x(10) = pressure loss between pressurant and propellant tanks, Pa
%       * x(11) = pressurant temperature, C
%       * x(12) = propellant reserve, %
%       * x(13) = trapped propellant, %
%       * x(14) = propellant loading uncertainty, %
%
% * ep: Environmental parameters
%       * ep(1) = delta V, m/s
%       * ep(2) = thrust, N
%       * ep(3) = tank maximum radius, m
%       * ep(4) = AOCS number of thrusters
%       * ep(5) = AOCS propellant mass, kg
%       * ep(6) = AOCS required thrust, N
%
%% Outputs:
% * M = total mass of the propulsion subsystem
% * P = electrical power required by the propulsion subsystem
% * info: 
%       * info.mass = mass inventory
%       * info.propellant = propellant characteristics
%
%% Author: Simone Alicino, 2012

M_sc = x(1);
Isp = x(2);
pc = x(3);
B = x(4);
epsilon = x(5);
if x(6) <= 1/3
    tank_material = 'aluminum';
elseif x(6) >= 2/3
    tank_material = 'titanium';
else
    tank_material = 'steel';
end
material_unc = 1 + x(7)/100;
if x(8) < 0.5
    pressurant = 'helium';
else
    pressurant = 'nitrogen';
end
p_pres0 = x(9);
p_loss = x(10);
T_pres = x(11);
prop_reserves = x(12)/100;
prop_trapped = x(13)/100;
prop_loading = x(14)/100;

Dv = ep(1);
F = ep(2);
r_tank_max = ep(3);
n_engines = ep(4);
aocs_propellant = ep(5);
F_aocs = ep(6);

g0 = 9.80665;
u = Isp*g0;
mdot = F/u;
stay_time = 0.01;
v_flow = 10;

% Select propulsion system and propellant from Isp
if (Isp >= 200) && (Isp <= 250)
    propellant.type = 'monopropellant';
    propellant.fuel = 'N2H4';
    propellant.density = 1.01e3;
    propellant.k = 1.27;
    propellant.M = 32.05;
    propellant.Tc = 1750;
elseif (Isp >= 250) && (Isp <= 335)
    propellant.type = 'bipropellant';
    propellant.ox = 'N2O4';
    propellant.fuel = 'N2H4';
    propellant.ox.density = 1.44e3;
    propellant.fuel.density = 1.01e3;
    propellant.MR = 1.1;
    propellant.density = (propellant.MR + 1)/(propellant.MR/propellant.ox.density + 1/propellant.fuel.density);
    propellant.k = 1.26;
    propellant.M = 1950;
    propellant.Tc = 3280;
elseif (Isp > 335) && (Isp <= 350)
    propellant.type = 'bipropellant';
    propellant.ox = 'LOX';
    propellant.fuel = 'RP-1';
    propellant.ox.density = 1.44e3;
    propellant.fuel.density = 0.91e3;
    propellant.MR = 2.3;
    propellant.density = (propellant.MR + 1)/(propellant.MR/propellant.ox.density + 1/propellant.fuel.density);
    propellant.k = 1.225;
    propellant.M = 2210;
    propellant.Tc = 3510;
elseif (Isp > 350) && (Isp <= 450)
    propellant.type = 'bipropellant';
    propellant.ox = 'LOX';
    propellant.fuel = 'LH2';
    propellant.ox.density = 1.142e3;
    propellant.fuel.density = 0.071e3;
    propellant.MR = 3.40;
    propellant.k = 1.26;
    propellant.M = 10;
    propellant.Tc = 2985;
end

% Propellant budget
M_propellant = M_sc*(1 - exp(-Dv/u)) + aocs_propellant;
% M_propellant = aocs_propellant;
M_reserves = prop_reserves*M_propellant;
M_usable = M_propellant + M_reserves;
M_trapped = prop_trapped*M_usable;
M_loadunc = prop_loading*M_usable;
M_totalprop = M_usable + M_trapped + M_loadunc;

ptank = pc + 0.2*pc + 0.5*propellant.density*v_flow^2 + p_loss;
switch propellant.type
    case 'monopropellant'
        [M_totaltank, V_ullage, V_usable] = tank(ptank, B, M_usable, prop_trapped, propellant.density, 'propellant', tank_material, material_unc, r_tank_max);
        M_pressurization = pressurization(ptank, p_pres0, p_loss, T_pres, V_ullage, V_usable, 'monopropellant', pressurant, r_tank_max, material_unc, B);
        M_mainengine = 0.4 + 0.0033*F;
        P_mainengine = interp1([1 2 4 22 275 440 3100],[13 14 15 30 35 45 168],F,'linear','extrap');
        M_rct = n_engines*(0.4 + 0.0033*F_aocs);
        P_rct = n_engines*interp1([1 2 4 22 275 440 3100],[13 14 15 30 35 45 168],F_aocs,'linear','extrap');
        
    case 'bipropellant'
        MR = propellant.MR;
        M_fuel = M_usable/(1 + MR);
        M_ox = M_usable - M_fuel;
        [M_tank.ox, V_ullage.ox, V_usable.ox] = tank(ptank, B, M_ox, prop_trapped, propellant.ox.density, 'propellant', tank_material, material_unc, r_tank_max);
        [M_tank.fuel, V_ullage.fuel, V_usable.fuel] = tank(ptank, B, M_fuel, prop_trapped, propellant.fuel.density, 'propellant', tank_material, material_unc, r_tank_max);
        M_totaltank = M_tank.ox + M_tank.fuel;
        M_pressurization = pressurization(ptank, p_pres0, p_loss, T_pres, V_ullage, V_usable, 'bipropellant', pressurant, r_tank_max, material_unc, B);
        M_mainengine = F/(0.0006098*F + 13.44)/g0;
        P_mainengine = interp1([22 110 490 890 4000],[5 36 45 46 70],F,'linear','extrap');
        M_rct = n_engines*(F_aocs/(0.0006098*F_aocs + 13.44)/g0);
        P_rct = n_engines*interp1([22 110 490 890 4000],[5 36 45 46 70],F_aocs,'linear','extrap');
end

% M_mainengine = thrust_chamber(pc, epsilon, F, mdot, stay_time, propellant, material_unc);
% M_rct = thrust_chamber(pc, epsilon, F_aocs, mdot, stay_time, propellant, material_unc);
% M_totalengine = M_mainengine.total + n_engines*M_rct.total;
M_totalengine = M_mainengine + M_rct;
P_totalengine = P_mainengine + P_rct;

% Valves, filters, ...
n_lines = strcmpi(propellant.type,'bipropellant') + 1;
m_filter = 0.1;
m_latch_valve = 0.4;
m_service_valve = 0.2;
m_orifice = 0.2;
m_transducer = 0.3;
m_pipes = 0.8;
M_feedsystem = n_lines*(m_pipes + m_filter + m_service_valve + m_orifice) + (2 + n_lines)*(m_latch_valve + m_transducer);
p_latch_valve = 0.11;
p_transducer = 0.28;
P_feedsystem = (2 + n_lines)*(p_latch_valve + p_transducer);

M = M_totalprop + M_totaltank + M_pressurization.total + M_totalengine + M_feedsystem;
P = P_totalengine + P_feedsystem;

mass.propellant.total = M_totalprop;
mass.propellant.required = M_propellant;
mass.propellant.reserves = M_reserves;
mass.propellant.usable = M_usable;
mass.tank = M_totaltank;
mass.pressurization = M_pressurization;
mass.mainengine = M_mainengine;
mass.rct = M_rct;
mass.feedsystem = M_feedsystem;
info.mass = mass;
info.propellant = propellant;

end

function [M_totaltank, V_ullage, V_usable] = tank(p, B, M_usable, prop_trapped, rho_propellant, tank_type, material, material_unc, r_tank_max)

% Material properties
aluminum.rho = 2800;
aluminum.sigma = 460e6;
steel.rho = 7860;
steel.sigma = 860e6;
titanium.rho = 4430;
titanium.sigma = 900e6;

switch lower(material)
    case 'aluminum'
        rho_tank = aluminum.rho;
        sigma = aluminum.sigma;
    case 'steel'
        rho_tank = steel.rho;
        sigma = steel.sigma;
    case 'titanium'
        rho_tank = titanium.rho;
        sigma = titanium.sigma;
end
rho_tank = rho_tank*material_unc;
sigma = sigma*material_unc;

V_usable = M_usable/rho_propellant;
V_propellant = (1+prop_trapped)*V_usable;

switch lower(tank_type)
    case 'propellant'
%         p_tank = 1.2*p;
        p_tank = p;
%         B = 3;
        V_ullage = V_usable/(B - 1);
        V_tank0 = V_propellant + V_ullage;
        % Diaphragm
        rho_diaphragm = 1519; % Teflon
        r_diaphragm = (0.75*V_tank0/pi)^(1/3);
        A_diaphragm = 2*pi*r_diaphragm^2;
        t_diaphragm = 0.002;
        V_diaphragm = A_diaphragm*t_diaphragm;
        M_diaphragm = V_diaphragm*rho_diaphragm;
        % Tank + diaphragm
        V_tank = V_tank0 + V_diaphragm;
        
    case 'pressurant'
        p_tank = p;
        V_ullage = 0;
        V_tank = V_propellant;
        M_diaphragm = 0;
end

% Spherical tank
cylinder = false;
r_tank = (0.75*V_tank/pi)^(1/3);
if r_tank > r_tank_max
    r_tank = r_tank_max;
    cylinder = true;
end
t_sphere = 0.5*p_tank*r_tank/sigma;
r_sphere = r_tank + t_sphere;
V_sphere = 4/3*pi*(r_sphere^3 - r_tank^3);
M_sphere = V_sphere*rho_tank;
M_tank = M_sphere + M_diaphragm;

% Cylindrical tank
if cylinder
    t_cylinder = p_tank*r_tank/sigma;
    r_cylinder = r_tank + t_cylinder;
    l_cylinder = (V_tank - 4/3*pi*r_tank^3)/(pi*r_tank^2);
    V_cylinder = pi*l_cylinder*(r_cylinder^2 - r_tank^2);
    M_cylinder = V_cylinder*rho_tank;
    M_tank = M_sphere + M_diaphragm + M_cylinder;
end

% Weldings and supports
M_girth = 2*pi*r_sphere*t_sphere*0.1*rho_tank;
M_penetrations = 2*pi*(0.075)^2*t_sphere*rho_tank;
M_welds = M_girth + M_penetrations;
M_supports = 0.02*(M_tank + M_welds);

M_totaltank = 2.5*(M_tank + M_welds + M_supports);

end

function [M_pres] = pressurization(p_tank, p0, p_loss, T_pres, V_ullage, V_usable, prop_type, pressurant, r_tank_max, material_unc, B)

% Pressurant properties
% T_pres = 20; % pressurant temperature, C
helium.R = 8314/4;
helium.T = 273.15 + T_pres;
nitrogen.R = 8314/28;
nitrogen.T = 273.15 + T_pres;

switch lower(pressurant)
    case 'helium'
        R = helium.R;
        T = helium.T;
    case 'nitrogen'
        R = nitrogen.R;
        T = nitrogen.T;
end

switch prop_type
    case 'monopropellant'
        % Blowdown
        M_pressurant = p_tank*V_ullage/(R*T);
        M_pres.total = M_pressurant;
        M_pres.pressurant = M_pressurant;
        
    case 'bipropellant'
        % Regulated Pressure
%         pf = p_tank + 7e5; % accounts for the regulator
        pf = p_tank;
        V_tankpres = p_tank*(V_usable.ox + V_usable.fuel)/(p0 - pf);
        M_tankpres = p0*V_tankpres/(R*T);
        M_ullagepres = p_tank*(V_ullage.ox + V_ullage.fuel)/(R*T);
        rho = p0/(R*T);
        M_tank = tank(p_tank, B, M_tankpres, 0, rho, 'pressurant', 'titanium', material_unc, r_tank_max);
        M_pres.total = M_tankpres + M_ullagepres + M_tank;
        M_pres.pressurant = M_tankpres + M_ullagepres;
        M_pres.tank = M_tank;
end

end

function [M_engine] = thrust_chamber(pc, epsilon, F, mdot, stay_time, propellant, material_unc)

k = propellant.k;
R = 8314/propellant.M;
Tc = propellant.Tc;

% Find exit pressure from area ratio and chamber pressure
Mach = @(M,e,k) -e+1/M*(2/(k+1)*(1+(k-1)*M^2/2))^((k+1)/2/(k-1));
Me = fzero(@(M) Mach(M,epsilon,k),3);
pe = pc*(1 + (k-1)/2*Me^2)^(k/(1-k));

Cf = sqrt(2*k^2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1 - (pe/pc)^((k-1)/k)) ) + (pe/pc)*epsilon;
cstar = sqrt(k*R*Tc)/(k*sqrt((2/(k+1))^((k+1)/(k-1))));

% Conical nozzle
At = F/(pc*Cf);
rt = sqrt(At/pi);
Ae = At*epsilon;
re = sqrt(Ae/pi);
Ln = (re - rt)/tand(15);

% Combustion chamber
steel.rho = 7860;
steel.sigma = 860e6;
carbon.rho = 1750;
carbon.sigma = 700e6;
inconel.rho = 8190;
inconel.sigma = 860e6;

rho = steel.rho*material_unc;
sigma = steel.sigma*material_unc;
time_residence = stay_time;
% mdot = pc*At/cstar;
Lstar = time_residence*R*Tc/cstar;
% ChTh = 2; % Sutton
% rc = ChTh*rt;
% Lc = Lstar/ChTh^2;
Mc = 0.2;
AcAt = 1/Mc*((2/(k+1))*(1+(k-1)/2*Mc^2))^((k+1)/(2*k-2));
rc = rt*sqrt(AcAt);
Lc = Lstar/AcAt;
tc = pc*rc/sigma;

Rc = rc + tc;
Rt = rt + tc;
Re = re + tc;

M_chamber = pi*Lc*rho*(Rc^2 - rc^2);
M_nozzle = 1/3*pi*Ln*carbon.rho*(Re^2 + Re*Rt + Rt^2 - re^2 - re*rt - rt^2);

% Injection plate
rho_p = propellant.density;
Ac = pi*rc^2;
% Vc = Ac*Lc;
% Cd = 0.99;
% Ai = rho_p*Vc/(time_residence*Cd*sqrt(0.4*pc*rho_p));
K = 1.5;
Dp = 0.2*pc;
Ai = mdot*sqrt(K/(2*rho_p*Dp));
% ti = 3*sqrt(4*Ai/pi); % injector length, L/D >= 3, Sutton
ti = sqrt(0.83*4.8*pc*rc/sigma); % thickness of a fully loaded flat circular plate, Roark
rho_inj = 1600;
M_injectors = rho_inj*ti*(Ac - Ai);
% M_injectors = 0;

M_engine.total = M_chamber + M_nozzle + M_injectors;
M_engine.chamber = M_chamber;
M_engine.nozzle = M_nozzle;
M_engine.injectors = M_injectors;

end
