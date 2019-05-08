function [M_batteries, battery] = batteries_sizing_new(E_req, cell_type, V_bus, par, eff)

% INPUT
% E_req: energy requirement.....(Wh)
% cell_type: descrete parameter
% V_bus.........................(V)
% eta: efficiency...............()

E_req = E_req/eff/par.fix.SoC_max;

% Eb = E_req/(eta.system.bdr*eta.system.lcl*eta.system.harn); % (Wh)
% DOD = Eb/Ecap;
% % Ecap = Eb/DOD;
% C = Ecap/V_bus;  %(Ah)


% capacities = [1.5 1.7 3.7 4.5];
% C1 = C;
% C1(C > capacities(end)) = C/2;
% capacity = sort([capacities C1]);
% selection = find(capacity == C1);
% selected = selection(end) - 1;
% selected(selected == 0) = 1;
% % selected = 1;

% cell.type = 'Li-ion';
% cell.efficiency = eta.system.batt;
% cell.voltage = 3.6;         % V
% cell.charge_voltage = 4.1;  % V

    
    % SSTL DATABASE

if cell_type <= .25
    % battery C
    cell.capacity = 1.5;         % Ah
    cell.energy = 5;             % Wh
    cell.spec_energy = 23.8095;  % Wh/kg
    cell.discharge_volt = 4.2;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.21;            % kg   
    cell.voltage = 4.2;          % V  
    cell.max_DoD = 0.75;         % Ah/Ah
    
elseif cell_type > .25 && cell_type <= .5
    % battery B
    cell.capacity = 1.7;         % Ah
    cell.energy = 5;             % Wh
    cell.spec_energy = 25;       % Wh/kg
    cell.discharge_volt = 4.2;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.20;            % kg   
    cell.voltage = 4.2;          % V    
    cell.max_DoD = 0.75;         % Ah/Ah
    
elseif cell_type > .5 && cell_type <= .75
    % battery D
    cell.capacity = 3.7;         % Ah
    cell.energy = 8.5;           % Wh
    cell.spec_energy = 36.9565;  % Wh/kg
    cell.discharge_volt = 4.1;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.23;            % kg   
    cell.voltage = 4.1;          % V 
    cell.max_DoD = 0.75;         % Ah/Ah
    
else    %if cell_type > .75 && cell_type <= 1
    % battery A
    cell.capacity = 4.5;         % Ah
    cell.energy = 16;            % Wh
    cell.spec_energy = 25.1429;  % Wh/kg
    cell.discharge_volt = 4.1;   % V
    cell.height = 250e-3;        % m (non sstl)
    cell.diameter = 53e-3;       % m (non sstl)
    cell.mass = 0.63;            % kg   
    cell.voltage = 4.1;          % V  
    cell.max_DoD = 0.8;         % Ah/Ah
end


battery.cells_series = ceil(V_bus/cell.voltage);
battery.cells_parallel = ceil(E_req/(cell.max_DoD*V_bus*cell.capacity))+1;
battery.cells_total = battery.cells_series*battery.cells_parallel;
battery.mass = battery.cells_total*cell.mass;
battery.cell = cell;

% battery.cell = cell;
% % battery.cells_series = ceil(V_bus/cell.discharge_volt);
% battery.required_capacity = C;
% battery.cells_series = ceil(V_bus/cell.voltage);
% battery.cells_parallel = ceil(C/cell.capacity);                               % + 1;   not failure at the moment
% battery.cells_total = battery.cells_series*battery.cells_parallel;
% battery.capacity = battery.cells_parallel*cell.capacity;
% battery.charge_voltage = battery.cells_series*cell.charge_voltage;
% battery.length = battery.cells_series*cell.diameter;
% battery.width = battery.cells_parallel*cell.diameter;
% battery.volume = battery.length*battery.width*cell.height;
% battery.mass = battery.cells_total*cell.mass;

M_batteries = battery.mass;


% battery.Ecap= Ecap;
end