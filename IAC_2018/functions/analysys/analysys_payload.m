% analysys of payload function

clear all; close all; clc

% NOMINAL DESIGN
tt              = 3000;
bit_depth       = 2;
camera_catalogs = 0.2;


% NOMINAL EPISTEMIC
h    = 600;    % altitude, epistemic (km)
eps0 = 2;      % elevation angle
delta_I = 0;   % delta inclination orbit (%)


% only change the inclination
ind_delta_I = 1;
for delta_I = 0:.001:1
    
    x = [tt bit_depth camera_catalogs];
    ep = [h eps0 delta_I];
    
    [M,P_imaging,P_idle,info] = space_payload_5subsystems(x,ep);
    plot_T_tt(ind_delta_I) = info.Tac_tot;
    
    ind_delta_I = ind_delta_I + 1;
end
plot([0:.001:1], plot_T_tt)




ind_delta_I = 1;
for delta_I = 0:.5:1
    
    ind_eps0 = 1;
    for eps0 = 0:2:4
        
        ind_h = 1;
        for h = 600:200:800
            
 
            ind_tt = 1;
            for tt = 0:6509
            
            x = [tt bit_depth camera_catalogs];
            ep = [h eps0 delta_I];
            
            [M,P_imaging,P_idle,info] = space_payload_5subsystems(x,ep);
            
            coverage(ind_delta_I, ind_eps0, ind_h) = info.coverage;
            Tac(ind_delta_I, ind_eps0, ind_h) = info.Tac_tot;
            
            plot_co(ind_h) = info.coverage;
            plot_T(ind_h) = info.Tac_tot;
            plot_T_tt(ind_tt) = info.Tac_tot;
%             save(strcat('payload_delta_I',num2str(delta_I),'_eps0',num2str(eps0),'_h',num2str(h)));

            ind_tt = ind_tt + 1;
            end
            hold on
            plot([0:6509], plot_T_tt)
            ind_h = ind_h + 1;
        end
        
        
        hold on            
%         plot([600:200:1200], plot_co)
%         plot([600:200:1200], plot_T)
        ind_eps0 = ind_eps0 + 1;
    end
    ind_delta_I = ind_delta_I + 1;
end




