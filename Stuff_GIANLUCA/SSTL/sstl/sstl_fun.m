function M_batteries = sstl_fun(d, u, par)

global nfevalglobal;
       nfevalglobal = nfevalglobal + 1;

%% fix parameters
day =  par.fix.day;
MJD0 = par.fix.MJD0;

%% uncertian parameters
eff = u(1);

%% design
V_bus = d(1);     %
cell_type = d(2); % [0:1]
%%

for ii=2:31
    u(ii-1) = u(ii);
end

ind = 1;
for i=1:5
  for j=1:6
    op(i,j)=u(ind);
    ind =ind+1;
  end
end

op(:,7) = [0; 24.7; 79.9; 114.7; 145.2];
op(:,8) = [0; .6; .8; .6; .3];


% op(:,1)
% op(:,2:end)


% %%%%%%%%%%%%%%%%%%%%%
% op =    [68500.3	0.902	22.81	86.63	180.1	0.0 	0.0     0.0;
%          73250.2	0.77866	9.12	86.79	180.06	180.08	24.7	0.6;
%          86065.5	0.513	1.09	85.96	180.81	180.84	79.9	0.8;
%          49646.4	0.15392	0.36	86.85	180.97	4.25	114.7	0.6;
%          42049      0.001	0.05	270     0       359.95	145.2	0.3];
% %%%%%%%%%%%%%%%%%%%%%     


% deg2rad
op(:,3:6) = deg2rad(op(:,3:6));

global unc
unc = [unc op];

t_on_off = battery_times(op,MJD0+day)
[e_req, e_sizing] = energy_required(t_on_off,eff);

% e_sizing

% mass = mask_batteries_sizing(d, e_sizing, eff);

global E
E = [E e_req];

% eta.system.bdr = eff;
% eta.system.lcl = eff;
% eta.system.harn = eff;
% eta.system.batt = eff;


% e_sizing.in(1) = E_cap*0.9;  % Battery start percentage: 90%
% 
% e_start(1) = e_sizing.in(1);
% for i=2:length(e_sizing.in)
%     e_start(i) = min(E_cap, e_start(i-1) - e_sizing.out(i-1) + e_sizing.in(i));
%     e_start(i) = max(0, e_start(i));
% end
% E_req = max(e_start - e_sizing.out);  % (Wh)


% [M_batteries, battery, DOD] = batteries_sizing(E_req, E_cap, V_bus, eta);
[M_batteries, battery] = batteries_sizing_new(e_req, cell_type, V_bus, par, eff);

return