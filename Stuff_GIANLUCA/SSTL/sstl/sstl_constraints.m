function [output_constraint, boh] = sstl_constraints(d, u, par)

% %% fix parameters
% day = par.fix.day;
% MJD0 = par.fix.MJD0;
% % ndays = 365;
% %% uncertian parameters
% eff = u(end);
% %% design
% V_bus = d(1);
% E_cap = d(2);
% %%
% 
% 
% ind = 1;
% for i=1:5
%   for j=1:6
%     op(i,j)=u(ind);
%     ind =ind+1;
%   end
% end
% 
% op(:,7) = [0; 24.7; 79.9; 114.7; 145.2];
% op(:,8) = [0; .6; .8; .6; .3];
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%
% % op =    [68500.3	0.902	22.81	86.63	180.1	0.0 	0.0     0.0;
% %          73250.2	0.77866	9.12	86.79	180.06	180.08	24.7	0.6;
% %          86065.5	0.513	1.09	85.96	180.81	180.84	79.9	0.8;
% %          49646.4	0.15392	0.36	86.85	180.97	4.25	114.7	0.6;
% %          42049      0.001	0.05	270     0       359.95	145.2	0.3];
% %%%%%%%%%%%%%%%%%%%%
% 
% 
% % deg2rad
% op(:,3:6) = deg2rad(op(:,3:6));
% 
% 
% t_on_off = battery_times(op,MJD0+day);
% [~, e_sizing] = energy_required(t_on_off,eff);
% 
% 
% 
% % mass = mask_batteries_sizing(d, e_sizing, eff);
% 
% 
% 
% eta.system.bdr = eff;
% eta.system.lcl = eff;
% eta.system.harn = eff;
% eta.system.batt = eff;
% 
% 
% 
% e_sizing.in(1) = E_cap;
% 
% e_start(1) = E_cap;
% for i=2:length(e_sizing.in)
%     e_start(i) = min(E_cap, e_start(i-1) - e_sizing.out(i-1) + e_sizing.in(i));
% end
% 
% E_req = max(e_start - e_sizing.out);
% [~, ~, DOD] = batteries_sizing(E_req, E_cap, V_bus, eta);
% 
% if DOD > 0.2 && DOD < 0.8
%     
%     output_constraint = 0;
% elseif DOD <= 0.2
%     
%     output_constraint = abs(DOD-0.2)*100;
% elseif DOD >= 0.8
%     
%     output_constraint = abs(DOD-0.8)*100;
% end


output_constraint = 0; 
boh = [];
return


