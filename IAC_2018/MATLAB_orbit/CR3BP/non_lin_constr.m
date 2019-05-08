function [c, ceq] = non_lin_constr(ctrl, pos0_vel0, t0, tf, posf_velf, Earth, Moon, L_star, T_star)

c = [];


[x, y, z, vx, vy, vz, t] = CR3BPIntegration3D_adim_imp_ctrl(pos0_vel0(1), pos0_vel0(2), pos0_vel0(3), ...
                                                            pos0_vel0(4), pos0_vel0(5), pos0_vel0(6), ...
                                                            ctrl, t0, tf, Earth, Moon, L_star, T_star);

 X_out = [x(end), y(end), z(end), vx(end), vy(end), vz(end)].';
 
 ceq = [];
 
 for index = 1:length(posf_velf)
     if isnan(posf_velf(index))
         continue;
     end
     
     ceq = [ceq; X_out(index)-posf_velf(index)];
 end