% function [M_batteries, battery] = batteries_sizing(Pecl, tecl, DOD, V_bus, eta)


function [M_batteries, battery, DOD] = batteries_sizing(E_req, Ecap, V_bus, eta)



% INPUT
% E_req: .......................(Wh)
% Ecap: capacity of the battery.(Wh)
% V_bus.........................(V)
% eta: efficiency...............()

% from Simone

Eb = E_req/(eta.system.bdr*eta.system.lcl*eta.system.harn); % (Wh)
DOD = Eb/Ecap;
% Ecap = Eb/DOD;
C = Ecap/V_bus;  %(Ah)

% capacities = [1.5 5.8 10 16 28 39 50];

capacities = [1.5 1.7 3.7 4.5];
C1 = C;
C1(C > capacities(end)) = C/2;
capacity = sort([capacities C1]);
selection = find(capacity == C1);
selected = selection(end) - 1;
selected(selected == 0) = 1;
% selected = 1;

cell.type = 'Li-ion';
cell.efficiency = eta.system.batt;
cell.voltage = 3.6;         % V
cell.charge_voltage = 4.1;  % V

% if selected <= 1
%     % Sony 18650 HC
%     cell.capacity = 1.5;        % Ah
%     cell.energy = 7;            % Wh
% %     cell.spec_energy = 146;     % Wh/kg
%     cell.spec_energy = 115;     % Wh/kg
%     cell.discharge_volt = 2.75; % V
%     cell.height = 65e-3;        % m
%     cell.diameter = 18.25e-3;   % m
%     cell.mass = 48e-3;          % kg
%     
% elseif selected == 2
%     % Saft MPS
%     cell.capacity = 5.8;        % Ah
%     cell.energy = 20;           % Wh
%     cell.spec_energy = 133;     % Wh/kg
%     cell.discharge_volt = 2.7;  % V
%     cell.height = 65e-3;        % m
%     cell.diameter = 18e-3;      % m
%     cell.mass = 0.15;           % kg
%     
% elseif selected == 3
%     % Saft VL 10E
%     cell.capacity = 10;         % Ah
%     cell.energy = 36;           % Wh
%     cell.spec_energy = 139;     % Wh/kg
%     cell.discharge_volt = 2.7;  % V
%     cell.height = 129e-3;       % m
%     cell.diameter = 34e-3;      % m
%     cell.mass = 0.25;           % kg
%     
% elseif selected == 4
%     % Saft VES 16
%     cell.capacity = 16;         % Ah
%     cell.energy = 16;           % Wh
%     cell.spec_energy = 155;     % Wh/kg
%     cell.discharge_volt = 2.57; % V
%     cell.height = 60e-3;        % m
%     cell.diameter = 33e-3;      % m
%     cell.mass = 0.155;          % kg
%     
% elseif selected == 5
%     % Saft VES 100
%     cell.capacity = 28;         % Ah
%     cell.energy = 100;          % Wh
%     cell.spec_energy = 118;     % Wh/kg
%     cell.discharge_volt = 2.57; % V
%     cell.height = 185e-3;       % m
%     cell.diameter = 54e-3;      % m
%     cell.mass = 0.81;           % kg
%     
% elseif selected == 6
%     % Saft VES 140
%     cell.capacity = 39;         % Ah
%     cell.energy = 140;          % Wh
%     cell.spec_energy = 126;     % Wh/kg
%     cell.discharge_volt = 2.7;  % V
%     cell.height = 250e-3;       % m
%     cell.diameter = 54e-3;      % m
%     cell.mass = 1.13;           % kg
%     
% elseif selected == 7
%     % Saft VES 180
%     cell.capacity = 50;         % Ah
%     cell.energy = 180;          % Wh
%     cell.spec_energy = 165;     % Wh/kg
%     cell.discharge_volt = 2.57; % V
%     cell.height = 250e-3;       % m
%     cell.diameter = 53e-3;      % m
%     cell.mass = 1.11;           % kg
    
    % SSTL DATABASE
    
if selected <= 1
    % battery C
    cell.capacity = 1.5;         % Ah
    cell.energy = 5;             % Wh
    cell.spec_energy = 23.8095;  % Wh/kg
    cell.discharge_volt = 4.2;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.21;            % kg   
    cell.voltage = 4.2;          % V    
    
elseif selected == 2
    % battery B
    cell.capacity = 1.7;         % Ah
    cell.energy = 5;             % Wh
    cell.spec_energy = 25;       % Wh/kg
    cell.discharge_volt = 4.2;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.20;            % kg   
    cell.voltage = 4.2;          % V    
    
elseif selected == 3
    % battery D
    cell.capacity = 3.7;         % Ah
    cell.energy = 8.5;           % Wh
    cell.spec_energy = 36.9565;  % Wh/kg
    cell.discharge_volt = 4.1;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.23;            % kg   
    cell.voltage = 4.1;          % V    
    
elseif selected == 4
    % battery A
    cell.capacity = 4.5;         % Ah
    cell.energy = 16;            % Wh
    cell.spec_energy = 25.1429;  % Wh/kg
    cell.discharge_volt = 4.1;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.63;            % kg   
    cell.voltage = 4.1;          % V   
end



battery.cell = cell;
% battery.cells_series = ceil(V_bus/cell.discharge_volt);
battery.required_capacity = C;
battery.cells_series = ceil(V_bus/cell.voltage);
battery.cells_parallel = ceil(C/cell.capacity);                               % + 1;   not failure at the moment
battery.cells_total = battery.cells_series*battery.cells_parallel;
battery.capacity = battery.cells_parallel*cell.capacity;
battery.charge_voltage = battery.cells_series*cell.charge_voltage;
battery.length = battery.cells_series*cell.diameter;
battery.width = battery.cells_parallel*cell.diameter;
battery.volume = battery.length*battery.width*cell.height;
battery.mass = battery.cells_total*cell.mass;

M_batteries = battery.mass;


battery.Ecap= Ecap;
end