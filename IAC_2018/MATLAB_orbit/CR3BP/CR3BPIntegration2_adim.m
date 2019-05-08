% Modified from Matlab code provided with the book Orbital Mechanics for
% Aerospace Engineers


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x_rot,y_rot,vx_rot,vy_rot,t] = CR3BPIntegration2_adim(x0_rot, y0_rot, vx0_rot, vy0_rot, time, Earth, Moon, L_star, T_star)
% ~~~~~~~~~~~~~~~~~~~
% Input: x0,y0,vx0 and vy0 in Earth-Moon rotating frame and tf, time of
% propagation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
  This program uses ode113 to solve the earth-moon restricted three-body problem
  Adapted from Matlab file from Curtis.

  days      - converts days to seconds
  Moon.R     - radius of the moon (km)
  Earth.R    - radius of the earth (km)
  r12       - distance from center of earth to center of moon (km)
  m1,m2     - masses of the earth and of the moon, respectively (kg)
  M         - total mass of the restricted 3-body system (kg)
  mu        - gravitational parameter of earth-moon system (km^3/s^2)
  Earth.mu,Moon.mu   - gravitational parameters of the earth and of the moon,
              respectively (km^3/s^2)
  pi_1,pi_2 - ratios of the earth mass and the moon mass, respectively,
              to the total earth-moon mass
  W         - angular velocity of moon around the earth (rad/s)
  x1,x2     - x-coordinates of the earth and of the moon, respectively,
              relative to the earth-moon barycenter (km)
  d0        - initial altitude of spacecraft (km)
  phi       - polar azimuth coordinate (degrees) of the spacecraft
              measured positive counterclockwise from the earth-moon line
  v0        - initial speed of spacecraft relative to rotating earth-moon
              system (km/s)
  gamma     - initial flight path angle (degrees)
  r0        - intial radial distance of spacecraft from the earth (km)
  x,y       - x and y coordinates of spacecraft in rotating earth-moon
              system (km)
  vx,vy     - x and y components of spacecraft velocity relative to
              rotating earth-moon system (km/s)
  f0        - column vector containing the initial valus of x, y, vx and vy
  t0,tf     - initial time and final times (s)
  t         - column vector of times at which the solution was computed
  f         - a matrix whose columns are:
              column 1: solution for x  at the times in t
              column 2: solution for y  at the times in t
              column 3: solution for vx at the times in t
              column 4: solution for vy at the times in t
  xf,yf     - x and y coordinates of spacecraft in rotating earth-moon
              system at tf
  vxf, vyf  - x and y components of spacecraft velocity relative to
              rotating earth-moon system at tf
  df        - distance from surface of the moon at tf
  vf        - relative speed at tf
  

  User M-functions required:  rkf45
  User subfunctions required: rates, circle
%}
% ---------------------------------------------

% clear all; close all; clc



% Initial position and velocity
f0     =  [x0_rot; y0_rot; vx0_rot; vy0_rot];


% Compute the trajectory:
[t,f]  = ode113(@(t,f)rates(t,f,Earth,Moon, L_star, T_star), time, f0);

x_rot      = f(:,1);
y_rot      = f(:,2);
vx_rot     = f(:,3);
vy_rot     = f(:,4);


return

% ~~~~~~~~~~~~~~~~~~~~~~~~
function dfdt = rates(t,f,Earth,Moon, L_star, T_star)
% ~~~~~~~~~~~~~~~~~~~~~~~~   
    %{
      This subfunction calculates the components of the relative acceleration
      for the restricted 3-body problem, using Equations 2.192a and 2.192b

      ax,ay - x and y components of relative acceleration (km/s^2)
      r1    - spacecraft distance from the earth (km)
      r2    - spacecraft  distance from the moon (km)
      f     - column vector containing x,  y,  vx and vy at time t
      dfdt  - column vector containing vx, vy, ax and ay at time t

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
    M      =  m1 + m2;
    
    % Dimensionless mass rations []
    pi1   =  m1 / M;
    pi2   =  m2 / M;
    
    mu = astro_constants(1) * M;
    mu = mu * T_star^2 / L_star^3;
    
    % Angular velocity
    AngVel  =  sqrt(mu / r12^3);
    
    % Position of mass 1 in the rotated reference frame
    x1     = -pi2 * r12;
    
    % Position of mass 2 in the rotated reference frame
    x2     =  pi1*r12;


    x      = f(1);
    y      = f(2);
    vx     = f(3);
    vy     = f(4);

    r1     = norm([x + pi2*r12, y]);
    r2     = norm([x - pi1*r12, y]);

    ax     =  2*AngVel*vy + AngVel^2*x - (Earth.mu * T_star^2 / L_star^3)*(x - x1)/r1^3 - ( Moon.mu * T_star^2 / L_star^3)*(x - x2)/r2^3;
    ay     = -2*AngVel*vx + AngVel^2*y - ((Earth.mu * T_star^2 / L_star^3)/r1^3 + ( Moon.mu * T_star^2 / L_star^3)/r2^3)*y;

    dfdt   = [vx; vy; ax; ay];
end 



end 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~