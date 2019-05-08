function [M,P,info] = space_power(x,ep)

%% space power: Electrical power system model
%
% [M,P,info] = space_power(x,ep)
%
%% Inputs:
% * x: Design parameters:
%       * x(1) = type of solar cell
%               Si = [0 0.33], GaAs28 = [0.33 0.66], GaAs30 = [0.66 1] 
%       * x(2) = required bus voltage, V [0,100]
%       * x(3) = eps configuration (det, mppt) [0,0.5) = det; [0.5,1] = mppt    
%              eps: electric power system, der: direct energy transfer, mppt: maximum point power tracking
%              (it is how the solar cells are connected to the battery, det=directly; mppt=using an interface that optimises the performance)
%       * x(4) = cell packing efficiency (or assembly factor) [0,1]
%       * x(5) = depth of discharge (it is how much the battery has been discharged) [0,1]
%       * x(6) = number of loads (to determine LCLs in PCDU) [1,50]
%              PCDU(Power Control and Distribution Unit):controls the power flow from the solar generator to the battery and distributes power, on command from the onboard computer, to the instruments, the heaters and the propulsion system
%              LCL(Latching Current Limiters)
%       * x(7) = harness mass factor (as a percentage of eps mass) (harness=cables and buses)
%       * x(8) = allowable voltage drop, %
%
% * ep: Environmental parameters:
%       * ep(1) = required daylight power, W
%       * ep(2) = required eclipse power, W
%       * ep(3) = orbital daylight time, h
%       * ep(4) = orbital eclipse time, h
%       * ep(5) = distance from Sun, AU
%       * ep(6) = worst case angle of incidence, deg
%       * ep(7) = uncertainty on array temperature, C
%       * ep(8) = satellite lifetime, yrs
%       * ep(9) = bus regulation (0: unregulated, 1: regulated)
%       * ep(10) = uncertainty on power requirements, %
%       * ep(11) = total radiation fluence
%       * ep(12) = s/c type (1: cubesat, 0: big one)
%
%% Outputs:
% * M = total mass of the electrical power subsystem
% * P = power produced by the electrical power subsystem
% * info:
%       * info.mass = mass inventory(structure with the mass of all the
%                      eps subsystems)
%       * info.array = array design(structure with all the info of the
%                      solar arrays)
%       * info.battery = battery design(structure with all the info of the
%                      batteries)
%       * info.modules = pcdu design(structure with all the info of the
%                      solar arrays)
%       * info.eta = efficiencies
%
%% Author: Simone Alicino, 2013
%  Reviewed by Aar�n del R�o Bellisco 2015

% Design variables
cell_type = x(1);
V_bus = x(2);
configuration = (x(3) >= 0.5); % 0 = det; 1 = mppt
eta.array.packing = x(4);
DOD = x(5);
n_loads = round(x(6));
k_harn = x(7);
v_drop = x(8); 

% Environmental parameters
Preq  = ep(1);
Pecl  = ep(2);
tsun  = ep(3); 
tecl  = ep(4);
range = ep(5);
theta = ep(6);
Tmargin  = ep(7);
life  = ep(8);
bus_type = logical(ep(9)>0.5); % 0 = unregulated; 1 = regulated
Pmargin = 1 + ep(10)/100;
MeV = ep(11);
Preq = Preq*Pmargin;
Pecl = Pecl*Pmargin;
cubesat = logical(ep(12)>0.5); % 0 = big; 1 = cubesat

% det
eta.system.sar = V_bus/(V_bus + 1.2);
% mppt
eta.system.sar(configuration) = interp1([20 30 40 50 70 100],[0.93 0.941 0.952 0.96 0.966 0.97],V_bus,'linear','extrap');
% unregulated
eta.system.bcr = 1; 
% regulated
eta.system.bcr(bus_type) = interp1([20 30 40 50 70 100],[0.90 0.92 0.935 0.95 0.96 0.97],V_bus,'linear','extrap');

eta.system.bdr = eta.system.bcr;

% Hardcoded parameters
eta.system.lcl = 0.99;
eta.system.batt = 0.96;
eta.system.harn = 1 - v_drop/100;

% Battery
[M_batteries, battery] = batteries(Pecl, tecl, DOD, V_bus, eta);

% Solar array
% V_sa = V_bus;
V_sa = battery.charge_voltage;
[M_array, P_array, array, eta] = solar_array(Preq, Pecl, tsun, tecl, range, theta, life, V_sa, cell_type, configuration, eta, MeV, Tmargin);

% PCDU
[M_pcdu, modules] = pcdu(configuration, bus_type, P_array, P_array-Preq, Pecl, n_loads, cubesat);

% Harness
% k_harn = density * resistivity * cable length^2 [kg-Ohm]
% 0.01 for 10m of cable
% 0.3  for 50m of cable
%  1   for 100m of cable
% M_harness = k_harn*P_array/V_bus^2/(1 - eta.system.harn);
% k_harn1 = k_harn/(1500/50^2/0.01);
% m_harn = k_harn1*P_array/V_bus^2/(1 - eta.system.harn);
m_harn = k_harn*0.02/(1 - eta.system.harn);
M_harness = m_harn*(M_array + M_pcdu + M_batteries);

% Total mass and power
M = M_array + M_pcdu + M_batteries + M_harness;
P = P_array;

mass.array = M_array;
mass.pcdu = M_pcdu;
mass.batteries = M_batteries;
mass.harness = M_harness;
mass.total = M;
info=struct;
info.mass = mass;
info.array = array;
info.battery = battery;
info.modules = modules;
info.eta = eta;

end

function [M_array, P_array, array, eta] = solar_array(Preq, Pecl, tsun, tecl, range, theta, life, V_bus, cell_type, configuration, eta, MeV, Tmargin)

fluences = [0 2.5e14 5e14 1e15];

if cell_type <= 0.33
        % Si  (Azurspace datasheet)
        cell.type = 'Si';
        cell.temp = 28;         % C
        cell.length = 74e-3;    % m
        cell.width = 32e-3;     % m
        cell.area = 23.61e-4;   % m^2
        cell.spec_mass = 32e-2; % kg/m^2
        cell.Voc = 0.628*interp1(fluences,[1,0.91,0.89,0.85],MeV,'linear','extrap');       % V
        cell.Isc = 458*cell.area*interp1(fluences,[1,0.88,0.85,0.76],MeV,'linear','extrap');% A
        cell.Vmp = 0.528*interp1(fluences,[1,0.91,0.89,0.84],MeV,'linear','extrap');       % V
        cell.Imp = 434*cell.area*interp1(fluences,[1,0.88,0.84,0.75],MeV,'linear','extrap');% A
        cell.efficiency = 0.169*interp1(fluences,[1,0.80,0.74,0.64],MeV,'linear','extrap'); % EOL
        cell.dVocdT = -1e-3*interp1(fluences,[2.02,2.14,2.17,2.20],MeV,'linear','extrap');    % V/C EOL
        cell.dIscdT = 10e-3*cell.area*interp1(fluences,[0.030,0.045,0.055,0.059],MeV,'linear','extrap');% A/C EOL
        cell.dVmpdT = -1e-3*interp1(fluences,[2.07,2.22,2.19,2.25],MeV,'linear','extrap');    % V/C EOL
        cell.dImpdT = 10e-3*cell.area*interp1(fluences,[0.004,0.023,0.024,0.035],MeV,'linear','extrap');% A/C EOL
        cell.absorptivity = 0.78;
        cell.emittance = 0.8;
        cell.degradation = 0.0375;% 1/year
        
        %---------------------------------
        % poisson distribution
        cell.interval_lambda = [0 20]; 
        % weibull distribution
        cell.interval_theta = [1800 2000];
        cell.interval_beta  = [0.2 0.6];
        %---------------------------------
        
elseif cell_type > 0.33 && cell_type <= 0.66
        % Triple Junction GaAs (Azurspace datasheet)
        cell.type = '3J GaAs';
        cell.temp = 28;         % C
        cell.length = 80e-3;    % m
        cell.width = 40e-3;     % m
        cell.area = 30.18e-4;   % m^2
        cell.spec_mass = 116e-2;% kg/m^2
        cell.Voc = 1e-3*interp1(fluences,[2662,2554,2528,2476],MeV,'linear','extrap');       % V
        cell.Isc = 1e-3*interp1(fluences,[505.4,500.9,500.8,485.8],MeV,'linear','extrap');% A
        cell.Vmp = 1e-3*interp1(fluences,[2365,2271,2226,2202],MeV,'linear','extrap');       % V
        cell.Imp = 1e-3*interp1(fluences,[487.0,482.1,472.4,457.8],MeV,'linear','extrap');% A
        cell.efficiency = interp1(fluences,[0.28,0.266,0.255,0.245],MeV,'linear','extrap'); % EOL
        cell.dVocdT = -1e-3*interp1(fluences,[6.0,6.4,6.2,6.3],MeV,'linear','extrap');    % V/C EOL
        cell.dIscdT = 1e-3*interp1(fluences,[0.32,0.33,0.31,0.39],MeV,'linear','extrap');% A/C EOL
        cell.dVmpdT = -1e-3*interp1(fluences,[6.1,6.8,6.3,6.4],MeV,'linear','extrap');    % V/C EOL
        cell.dImpdT = 1e-3*interp1(fluences,[0.28,0.36,0.20,0.29],MeV,'linear','extrap');% A/C EOL
        cell.absorptivity = 0.91;
        cell.emittance = 0.84;
        cell.degradation = 0.0275;% 1/year
        %---------------------------------
        % poisson distribution
        cell.interval_lambda = [10 30];
        % weibull distribution
        cell.interval_theta = [1900 2100];
        cell.interval_beta  = [0.25 0.65];
        %---------------------------------
        
elseif cell_type > 0.66
        % Triple Junction GaAs (Azurspace datasheet)
        cell.type = '3J GaAs';
        cell.temp = 28;         % C
        cell.length = 80e-3;    % m
        cell.width = 40e-3;     % m
        cell.area = 30.18e-4;   % m^2
        cell.spec_mass = 86e-2;% kg/m^2
        cell.Voc = 1e-3*interp1(fluences,[2690,2560,2514,2468],MeV,'linear','extrap');       % V
        cell.Isc = 1e-3*interp1(fluences,[519.6,517.1,514.3,501.3],MeV,'linear','extrap');% A
        cell.Vmp = 1e-3*interp1(fluences,[2409,2277,2229,2191],MeV,'linear','extrap');       % V
        cell.Imp = 1e-3*interp1(fluences,[502.9,499.2,493.4,477.6],MeV,'linear','extrap');% A
        cell.efficiency = interp1(fluences,[0.293,0.276,0.267,0.254],MeV,'linear','extrap'); % EOL
        cell.dVocdT = -1e-3*interp1(fluences,[6.0,6.1,6.2,6.3],MeV,'linear','extrap');    % V/C EOL
        cell.dIscdT = 1e-3*interp1(fluences,[0.32,0.35,0.31,0.39],MeV,'linear','extrap');% A/C EOL
        cell.dVmpdT = -1e-3*interp1(fluences,[6.1,6.2,6.3,6.4],MeV,'linear','extrap');    % V/C EOL
        cell.dImpdT = 1e-3*interp1(fluences,[0.28,0.27,0.20,0.29],MeV,'linear','extrap');% A/C EOL
        cell.absorptivity = 0.91;
        cell.emittance = 0.84;
        cell.degradation = 0.0275;% 1/year
        %---------------------------------
        % poisson distribution
        cell.interval_lambda = [20 40];
        % weibull distribution
        cell.interval_theta = [1950 2191];
        cell.interval_beta  = [0.4 0.45];
        %---------------------------------
        
end
substrate_spec_mass = 3.2; % kg/m^2

% Efficiencies and losses
eta.array.pointing = kellycosd(theta);
eta.array.degr = (1 - cell.degradation)^life;
eta.array.range = 1/range^2;
% eta.array.otherfactors = 0.98*0.99*0.99;
eta.array.otherfactors = 0.99*0.97*0.99*0.9875*0.99;
eta.array.total = eta.array.packing*eta.array.pointing*eta.array.degr*eta.array.range*eta.array.otherfactors;

Xd = eta.system.lcl*eta.system.harn;
Xe = eta.system.bcr*eta.system.bdr*eta.system.batt*eta.system.lcl*eta.system.harn;

P_array = (((Pecl * tecl)/Xe) + ((Preq * tsun)/Xd)) / tsun;
V_array = V_bus/eta.system.sar;

% Solar array temperature (Agrawal, p. 285)
G = 1365;
sigma = 5.67e-8;
epsilon_B = 0.7;
alpha_OSR = 0.2;
P_cell = G*cell.efficiency*eta.array.total;
A_est = P_array/P_cell;
alpha_SE = (cell.absorptivity*eta.array.packing + alpha_OSR*(1-eta.array.packing)) - cell.efficiency;
S = G*eta.array.range*eta.array.pointing;
Tmax = nthroot( (alpha_SE*A_est*S)/((cell.emittance + epsilon_B)*A_est*sigma),4) - 273.15 + Tmargin;

% Operative voltage and current
Vmp = cell.Vmp + cell.dVmpdT*(Tmax - cell.temp);
Voc = cell.Voc + cell.dVocdT*(Tmax - cell.temp);
Imp = cell.Imp + cell.dImpdT*(Tmax - cell.temp);
Isc = cell.Isc + cell.dIscdT*(Tmax - cell.temp);
Ca = (Vmp/Voc - 1)/log(1 - Imp/Isc);
Cb = (1 - Imp/Isc)*exp(-Vmp/Voc/Ca);
cell.Vop = 0.97*Vmp;
% cell.Iop = Isc*( 1 - Cb*(exp(-cell.Vop/Voc/Ca) - 1) );

% Array sizing
V_cell = cell.Vop;
V_cell(logical(configuration)) = Vmp;
I_cell = Isc*( 1 - Cb*(exp(-V_cell/Voc/Ca) - 1) );
P_cell = I_cell*V_cell*eta.array.total;
Nc = ceil(P_array / P_cell);


%--------------------------------------------------------------------------
% poisson distribution
array.interval_lambda = cell.interval_lambda;
% weibull distribution
array.interval_theta = cell.interval_theta;
array.interval_beta = cell.interval_beta;
%--------------------------------------------------------------------------

array.cell = cell;
array.cells_series = ceil(V_array / V_cell);        % Number of cells in series
array.cells_parallel = ceil(Nc / array.cells_series) + 1;  % Number of strings (cells in parallel) + 1 redundant
array.cells_total = array.cells_series*array.cells_parallel;

array.length = array.cells_parallel*cell.length;
array.width = array.cells_series*cell.width;
array.area = cell.area*array.cells_total;
array.mass = array.area*(cell.spec_mass + substrate_spec_mass);
array.voltage = array.cells_series*V_cell;
array.current = array.cells_parallel*I_cell*eta.array.total;
array.power = array.voltage*array.current;
array.temperature = Tmax;

M_array = array.mass;
P_array = array.power;

end

function [M_pcdu, modules] = pcdu(configuration, bus_type, Parray, Pcharge, Pdischarge, n_loads, cubesat)

% m_module = 0.5; % average mass of PCDU module, kg
m_module = 0.75; % average mass of PCDU module, kg

Psar = 500;
Pbcr = 300;
Pbdr = 600;
Pmppt = Psar;
Plcl = 900;
Pheaters = 400;

if cubesat
        modules.ctrl = (Parray/500)*m_module;
        modules.tmtc = (Parray/500)*m_module;
        modules.bcr = bus_type*(Pcharge/Pbcr)*m_module;
        modules.bdr = bus_type*(Pdischarge/Pbdr)*m_module;
        modules.mppt = configuration*(Parray/Pmppt)*m_module;%*ceil(Parray/Pmppt);
        modules.sar = (Parray/Psar)*m_module;
        % modules.pyro = ceil(n_plines/25);
        % modules.heaters = ceil(n_hlines/16);
        modules.lcl = ceil(n_loads/16)*(Parray/Plcl)*m_module;
        casing = 1;
else
        modules.ctrl = 1;
        modules.tmtc = 1;
        modules.bcr = bus_type*ceil(Pcharge/Pbcr);
        modules.bdr = bus_type*(ceil(Pdischarge/Pbdr) + 1);
        modules.mppt = configuration*1;%*ceil(Parray/Pmppt);
        modules.sar = ceil(Parray/Psar) + 0*1;
        % modules.pyro = ceil(n_plines/25);
        % modules.heaters = ceil(n_hlines/16);
        modules.pyro = 2;
        modules.heaters = 2;
        modules.lcl = ceil(n_loads/16);%*ceil(Parray/Plcl);
        casing = m_module + 0.27;
end

modules.total = sum(cell2mat(struct2cell(modules)));
modules.mass = modules.total*casing;

M_pcdu = modules.mass;

%pcdu dimensions
if cubesat
%(based on the Cubesat pcdu of Clyde Space Ltd)
%http://www.clyde-space.com/cubesat_shop/power_distribution_and_protection/95_cubesat-power-distribution-module
    modules.dimensions.length=0.09;
    modules.dimensions.width=0.09;
    modules.dimensions.height=0.024;
else
%(based on the LEO pcdu of Surrey Satellite Technology Ltd
% http://www.sstl.co.uk/getattachment/fb984708-a8d0-4cc5-9253-49541bf85151/Power-Control-Distribution-Unit-PCDU)
    modules.dimensions.lenght=0.335;
    modules.dimensions.width=0.305;
    longitudinal_density=4.4/0.078;
    modules.dimensions.height=modules.mass/longitudinal_density;
end



end

function [M_batteries, battery] = batteries(Pecl, tecl, DOD, V_bus, eta)


Eb = Pecl*tecl/(eta.system.bdr*eta.system.lcl*eta.system.harn);
Ecap = Eb/DOD;
C = Ecap/V_bus;

capacities = [1.5 5.8 10 16 28 39 50];
C1 = C;
C1(C > capacities(end)) = C/2;
capacity = sort([capacities C1]);
% selection = find(capacity == C1);
% selected = selection(end) - 1;
% selected(selected == 0) = 1;
selected = 1;

cell.type = 'Li-ion';
cell.efficiency = eta.system.batt;
cell.voltage = 3.6;         % V
cell.charge_voltage = 4.1;  % V

if selected <= 1
    % Sony 18650 HC
    cell.capacity = 1.5;        % Ah
    cell.energy = 7;            % Wh
%     cell.spec_energy = 146;     % Wh/kg
    cell.spec_energy = 115;     % Wh/kg
    cell.discharge_volt = 2.75; % V
    cell.height = 65e-3;        % m
    cell.diameter = 18.25e-3;   % m
    cell.mass = 48e-3;          % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [2 4];
    % weibull distribution
    cell.interval_theta = [7000 8000];
    cell.interval_beta  = [0.6 0.8];
    %---------------------------------
    
elseif selected == 2
    % Saft MPS
    cell.capacity = 5.8;        % Ah
    cell.energy = 20;           % Wh
    cell.spec_energy = 133;     % Wh/kg
    cell.discharge_volt = 2.7;  % V
    cell.height = 65e-3;        % m
    cell.diameter = 18e-3;      % m
    cell.mass = 0.15;           % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [20 40];
    % weibull distribution
    cell.interval_theta = [7500 8500];
    cell.interval_beta  = [0.59 0.9];
    %---------------------------------
        
elseif selected == 3
    % Saft VL 10E
    cell.capacity = 10;         % Ah
    cell.energy = 36;           % Wh
    cell.spec_energy = 139;     % Wh/kg
    cell.discharge_volt = 2.7;  % V
    cell.height = 129e-3;       % m
    cell.diameter = 34e-3;      % m
    cell.mass = 0.25;           % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [10 15];
    % weibull distribution
    cell.interval_theta = [6500 7600];
    cell.interval_beta  = [0.61 0.88];
    %---------------------------------
    
elseif selected == 4
    % Saft VES 16
    cell.capacity = 16;         % Ah
    cell.energy = 16;           % Wh
    cell.spec_energy = 155;     % Wh/kg
    cell.discharge_volt = 2.57; % V
    cell.height = 60e-3;        % m
    cell.diameter = 33e-3;      % m
    cell.mass = 0.155;          % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [20 21];
    % weibull distribution
    cell.interval_theta = [7400 8100];
    cell.interval_beta  = [0.62 0.85];
    %---------------------------------
    
elseif selected == 5
    % Saft VES 100
    cell.capacity = 28;         % Ah
    cell.energy = 100;          % Wh
    cell.spec_energy = 118;     % Wh/kg
    cell.discharge_volt = 2.57; % V
    cell.height = 185e-3;       % m
    cell.diameter = 54e-3;      % m
    cell.mass = 0.81;           % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [6 58];
    % weibull distribution
    cell.interval_theta = [5000 10000];
    cell.interval_beta  = [0.66 0.89];
    %---------------------------------
    
elseif selected == 6
    % Saft VES 140
    cell.capacity = 39;         % Ah
    cell.energy = 140;          % Wh
    cell.spec_energy = 126;     % Wh/kg
    cell.discharge_volt = 2.7;  % V
    cell.height = 250e-3;       % m
    cell.diameter = 54e-3;      % m
    cell.mass = 1.13;           % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [32 44];
    % weibull distribution
    cell.interval_theta = [7000 8057];
    cell.interval_beta  = [0.64 0.87];
    %---------------------------------
    
elseif selected == 7
    % Saft VES 180
    cell.capacity = 50;         % Ah
    cell.energy = 180;          % Wh
    cell.spec_energy = 165;     % Wh/kg
    cell.discharge_volt = 2.57; % V
    cell.height = 250e-3;       % m
    cell.diameter = 53e-3;      % m
    cell.mass = 1.11;           % kg
    %---------------------------------
    % poisson distribution
    cell.interval_lambda = [0 10];
    % weibull distribution
    cell.interval_theta = [7000 7733];
    cell.interval_beta  = [0.7 0.84];
    %---------------------------------
    
end


%--------------------------------------------------------------------------
% poisson distribution
battery.interval_lambda = cell.interval_lambda;
% weibull distribution
battery.interval_theta = cell.interval_theta;
battery.interval_beta = cell.interval_beta;
%--------------------------------------------------------------------------
    
battery.cell = cell;
% battery.cells_series = ceil(V_bus/cell.discharge_volt);
battery.required_capacity = C;
battery.cells_series = ceil(V_bus/cell.voltage);
battery.cells_parallel = ceil(C/cell.capacity) + 1;
battery.cells_total = battery.cells_series*battery.cells_parallel;
battery.capacity = battery.cells_parallel*cell.capacity;
battery.charge_voltage = battery.cells_series*cell.charge_voltage;
battery.length = battery.cells_series*cell.diameter;
battery.width = battery.cells_parallel*cell.diameter;
battery.volume = battery.length*battery.width*cell.height;
battery.mass = battery.cells_total*cell.mass;

M_batteries = battery.mass;

end

function  y = kellycosd(x)

% Kelly cosine, x in degrees

y = cosd(x);
y(x >= 50) = interp1([50 60 80 85],[.635 .450 .1 0],x,'linear',0);

end
