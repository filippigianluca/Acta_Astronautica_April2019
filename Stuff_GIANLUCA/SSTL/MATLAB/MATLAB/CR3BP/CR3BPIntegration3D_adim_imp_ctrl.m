% Modified from Matlab code provided with the book Orbital Mechanics for
% Aerospace Engineers


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x_rot,y_rot,z_rot,vx_rot,vy_rot,vz_rot,tf_out] = CR3BPIntegration3D_adim_imp_ctrl(x0_rot, y0_rot, z0_rot,...
                                                                                        vx0_rot, vy0_rot, vz0_rot,...
                                                                                        ctrl,...
                                                                                        t0, tf, Earth, Moon, L_star, T_star)
% ~~~~~~~~~~~~~~~~~~~
% Input: x0,y0,z0,vx0,vy0 and vz0 in Earth-Moon rotating frame, ctrl
%        and tf, time of propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
  This program uses ode113 to solve the earth-moon restricted three-body problem
  Adapted from Matlab file from Curtis.

  days        - converts days to seconds
  Moon.R      - radius of the moon (km)
  Earth.R     - radius of the earth (km)
  r12         - distance from center of earth to center of moon (km)
  m1,m2       - masses of the earth and of the moon, respectively (kg)
  M           - total mass of the restricted 3-body system (kg)
  mu          - gravitational parameter of earth-moon system (km^3/s^2)
  Earth.mu    - gravitational parameters of the earth (km^3/s^2)
  Moon.mu     - gravitational parameters of the moon (km^3/s^2)
  pi_1,pi_2   - ratios of the earth mass and the moon mass, respectively,
                to the total earth-moon mass
  W           - angular velocity of moon around the earth (rad/s)
  x1,x2       - x-coordinates of the earth and of the moon, respectively,
                relative to the earth-moon barycenter (km)
  d0          - initial altitude of spacecraft (km)
  phi         - polar azimuth coordinate (degrees) of the spacecraft
                measured positive counterclockwise from the earth-moon line
  v0          - initial speed of spacecraft relative to rotating earth-moon
                system (km/s)
  gamma       - initial flight path angle (degrees)
  r0          - intial radial distance of spacecraft from the earth (km)
  x,y,z       - x, y and z coordinates of spacecraft in rotating earth-moon
                system (km)
  vx,vy,vz    - x, y and z components of spacecraft velocity relative to
                rotating earth-moon system (km/s)
  f0          - column vector containing the initial values of x, y, z,
                vx, vy and vz
  t0,tf       - initial time and final times (s)
  t           - column vector of times at which the solution was computed
  f           - a matrix whose columns are:
                column 1: solution for x  at the times in t
                column 2: solution for y  at the times in t
                column 3: solution for z  at the times in t
                column 4: solution for vx at the times in t
                column 5: solution for vy at the times in t
                column 6: solution for vz at the times in t
  xf,yf,zf    - x, y and z coordinates of spacecraft in rotating earth-moon
                system at tf
  vxf,vyf,vzf - x, y and z components of spacecraft velocity relative to
                rotating earth-moon system at tf
  df          - distance from surface of the moon at tf
  vf          - relative speed at tf
  

  User M-functions required:  rkf45
  User subfunctions required: rates, circle
%}
% ---------------------------------------------

% clear all; close all; clc



% Initial position and velocity
f0 = [x0_rot; y0_rot; z0_rot; vx0_rot; vy0_rot; vz0_rot];


% Control extraction
if isempty(tf)
    time = [t0 ctrl(end)];
    ctrl = ctrl(1:end-1);
else
    time = [t0 tf];
end

ctrl_seq = reshape(ctrl, 4, []).';
t_ctrl = ctrl_seq(:,1);


% Integration options
options = odeset('RelTol',1e-10,'AbsTol',1e-10);


% Compute the trajectory:
t_int = [time(1) t_ctrl(1)];
f0_int = f0;

tf_out = time(1);
ff = f0.';

for index = 1:length(t_ctrl)+1
    
    [tf_int, ff_int] = ode113(@(t,f)rates_3D_imp_ctrl(t, f, Earth, Moon, L_star, T_star), t_int, f0_int, options);
    
    tf_out = [tf_out; tf_int(2:end)];
    ff = [ff; ff_int(2:end,:)];
    
    if index == length(t_ctrl)
        t_int = [tf_out(end) time(end)];
    elseif index == length(t_ctrl)+1
        break;
    else
        t_int = [tf_out(end) t_ctrl(index+1)];
        ff(end, 4:6) = ff(end, 4:6)+ctrl_seq(index, 2:4);
    end
    
    f0_int = ff(end,:).';
    
end

x_rot  = ff(:,1);
y_rot  = ff(:,2);
z_rot  = ff(:,3);
vx_rot = ff(:,4);
vy_rot = ff(:,5);
vz_rot = ff(:,6);


return

% ~~~~~~~~~~~~~~~~~~~~~~~~
function dfdt = rates_3D_imp_ctrl(t, f, Earth, Moon, L_star, T_star)
% ~~~~~~~~~~~~~~~~~~~~~~~~   
    %{
      This subfunction calculates the components of the relative acceleration
      for the restricted 3-body problem, using Equations 2.192a and 2.192b

      ax,ay,az - x, y and z components of relative acceleration (km/s^2)
      r1    - spacecraft distance from the earth (km)
      r2    - spacecraft  distance from the moon (km)
      f     - column vector containing x, y, z, vx, vy and vz at time t
      dfdt  - column vector containing vx, vy, vz, ax, ay and az at time t

      All other variables are defined above.

      User M-functions required: none
    %}
    % ------------------------    
    % Earth-Moon distance [km] and [L_star]
    r12    =  Moon.mean_distance;
    r12 = r12 / L_star;
    
    % Earth and Moon masses [kg]
    m1 = Earth.mass;
    m2 = Moon.mass;
    
    % Total mass of the system [kg]
    M = m1 + m2;
    
    % Dimensionless mass rations []
    pi1 = m1 / M;
    pi2 = m2 / M;
    
    mu = astro_constants(1) * M;
    mu = mu * T_star^2 / L_star^3;
    
    % Angular velocity
    AngVel = sqrt(mu / r12^3);
    
    % Position of mass 1 in the rotated reference frame
    x1 = -pi2 * r12;
    
    % Position of mass 2 in the rotated reference frame
    x2 = pi1*r12;


    x    = f(1);
    y    = f(2);
    z    = f(3);
    vx   = f(4);
    vy   = f(5);
    vz   = f(6);

    r1   = norm([x + pi2*r12, y, z]);
    r2   = norm([x - pi1*r12, y, z]);

    ax   =  2*AngVel*vy + AngVel^2*x - (Earth.mu * T_star^2 / L_star^3)*(x - x1)/r1^3 - ( Moon.mu * T_star^2 / L_star^3)*(x - x2)/r2^3;
    ay   = -2*AngVel*vx + AngVel^2*y - ((Earth.mu * T_star^2 / L_star^3)/r1^3 + ( Moon.mu * T_star^2 / L_star^3)/r2^3)*y;
    az   = -((Earth.mu * T_star^2 / L_star^3)/r1^3 + ( Moon.mu * T_star^2 / L_star^3)/r2^3)*z;

    dfdt = [vx; vy; vz; ax; ay; az];
end 

end 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~