
%% cubesat

kepel = [7000000, 0.01, 98*deg2rad, 0, 35*deg2rad, 0];   

% Orbit state vector:
stat = kepel_statvec(kepel);
    
% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel); % just to compare with the numeric orbit propagator



% Propagation time in seconds:
tstart  = 0;        % initial time (sec)
tstep   = 1;        % step time (sec)
tend    = 6000;     % end time (~1 orbit)
n       = fix(tend/tstep);  % data size




%% sun
kepel_sun = [7000000000, 0, 0, 0, 0, 0];   
% Orbit state vector:
stat_sun = kepel_statvec(kepe_sun);
% Compute the variations in keplerian elements due to the Earth oblateness
delk_sun = delkep(kepe_sun); % just to compare with the numeric orbit propagator







% Orbit propagation
for t = tstart:tstep:tend

    % Analytical orbit propagation
    kp_an   = kepel + delk*t;
    
    % Convert from keplerian elements to state vector
    sv_an   = kepel_statvec(kp_an);
    xi_an   = sv_an(1:3)';
    vi_an   = sv_an(4:6)';

    % Orbit reference frame rotation matrix
    c_i_o   = orbital_to_inertial_matrix(kp_an);
    
    
    r_xo(:, ic) = c_i_o'*(xi_nm - xi_an)/1000;
    r_vo(:, ic) = c_i_o'*(vi_nm - vi_an);
     
    
  %% sun
     % Analytical orbit propagation
    kp_an_sun   = kepel_sun + delk_sun*t;
    
    % Convert from keplerian elements to state vector
    sv_an_sun   = kepel_statvec(kp_an_sun);
    xi_an_sun   = sv_an_sun(1:3)';
    vi_an_sun   = sv_an_sun(4:6)';

    % Orbit reference frame rotation matrix
    c_i_o_sun   = orbital_to_inertial_matrix(kp_an_sun);
    
    
    r_xo_sun(:, ic) = c_i_o'*(xi_nm - xi_an)/1000;
    r_vo_sun(:, ic) = c_i_o'*(vi_nm - vi_an); 
  
    
    
    
end



