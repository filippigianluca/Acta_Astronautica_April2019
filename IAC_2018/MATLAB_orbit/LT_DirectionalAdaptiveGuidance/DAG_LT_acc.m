function f_LT = DAG_LT_acc(x, DAG_parameters, constants)

% Input: x ->
%        DAG_parameters ->
%        constants -> 

% Output: f_LT -> low-thrust acceleration vector

% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

x0          = DAG_parameters.x0;
x_target    = DAG_parameters.x_target;
target_flag = DAG_parameters.target_flag;
tol         = DAG_parameters.tol;

% Current orbital elements
a    = x(1);
e    = x(2);
incl = x(3);
RAAN = x(4);
omega = x(5);
theta = x(6);

r = a * (1-e^2) / (1 + e * cos(theta));
v = sqrt(2 * constants.mu / r - constants.mu / a);

% Initial orbital elements
a0    = x0(1);
e0    = x0(2);
incl0 = x0(3);
RAAN0 = x0(4);
omega0 = x0(5);
theta0 = x0(6);

% Target orbital elements
a_target    = x_target(1);
e_target    = x_target(2);
incl_target = x_target(3);
RAAN_target = x_target(4);
omega_target = x_target(5);
theta_target = x_target(6);

cos_E = (e + cos(theta)) / (1 + e * cos(theta));

%% Semimajor axis

    
alpha_a = atan2(e * cos(theta), 1 + e * cos(theta));
beta_a  = 0;

f_a = [cos(beta_a) * sin(alpha_a); ...
    cos(beta_a) * cos(alpha_a); ...
    sin(beta_a)];

if target_flag(1)    
    
    R_a = ( a_target - a ) / abs(a_target - a0);
    if abs(R_a) < 1e-4
        R_a = 0;
    end
    
    eta_a = v * sqrt(a * (1-e) / (constants.mu * (1 + e)));
else
    R_a = 0;
    eta_a = 0;
end

%% Eccentricity

alpha_e = atan2(sin(theta), cos(theta) + cos_E);
beta_e  = 0;

f_e = [cos(beta_e) * sin(alpha_e); ...
       cos(beta_e) * cos(alpha_e); ...
       sin(beta_e)];
   
if target_flag(2)
   R_e = ( e_target - e ) / abs(e_target - e0);
   if abs(R_e) < 1e-6
       R_e = 0;
   end

   eta_e = (1 + 2 * e * cos(theta) + cos(theta)^2) / (2 * (1 + e * cos(theta)));

else
   R_e = 0;
   eta_e = 0;
end

   
%% Inclinaton 

alpha_incl = 0;
beta_incl  = sign( cos(omega+theta) ) * pi/2;

f_incl = [cos(beta_incl) * sin(alpha_incl); ...
    cos(beta_incl) * cos(alpha_incl); ...
    sin(beta_incl)];

if target_flag(3)
    R_incl = ( incl_target - incl ) / abs(incl_target - incl0);
    if abs(R_incl) < 1e-3
        R_incl = 0;
    end
    
    eta_incl = abs( cos( omega+theta) ) / (1 + e * cos(theta)) * ...
        (sqrt(1 - e^2 * sin(omega)^2) - e * abs( cos(omega) ) );
else
    R_incl   = 0;
    eta_incl = 0;
end


%% Right ascension

alpha_RAAN = 0;
beta_RAAN  = sign( sin(omega+theta) ) * pi/2;

f_RAAN = [cos(beta_RAAN) * sin(alpha_RAAN); ...
       cos(beta_RAAN) * cos(alpha_RAAN); ...
       sin(beta_RAAN)];

if target_flag(4)
   R_RAAN = ( RAAN_target - RAAN ) / abs(RAAN_target - RAAN0);
   if abs(R_RAAN) < 1e-4
       R_RAAN = 0;
   end

   eta_RAAN = abs( sin( omega+theta) ) / (1 + e * cos(theta)) * ...
       (sqrt(1 - e^2 * sin(omega)^2) - e * abs( sin(omega) ) );
else
   R_RAAN = 0;
   eta_RAAN = 0;
end


%% Argument of the perigee
alpha_omega = atan2( (1 + e * cos(theta)) * cos(theta), ...
                      ((2 + e * cos(theta)) * sin(theta)) );
     
beta_omega = atan2( e * cot(incl) * sin( omega+theta ) , ...
    ( sin(alpha_omega - theta) * (1 + e * cos(theta)) - cos(alpha_omega) * sin(theta)) );

f_omega = [cos(beta_omega) * sin(alpha_omega); ...
       cos(beta_omega) * cos(alpha_omega); ...
       sin(beta_omega)];
   
   
if target_flag(5)
   R_omega = ( omega_target - omega ) / abs(omega_target - omega0);
   if abs(R_omega) < 1e-4
       R_omega = 0;
   end

   theta_omega_max = acos( ( (1-e^2) / (2 * e^3) + ...
       sqrt(0.25 * (1-e^2) / (2 * e^3) +1/27) )^(1/3) + ...
       ( -(1-e^2) / (2 * e^3) + ...
       sqrt(0.25 * (1-e^2) / (2 * e^3) +1/27) )^(1/3) - 1 /e );

   eta_omega = (1 + sin(theta)^2) / (1 + e * cos(theta)) * ...
       (1 + e * cos(theta_omega_max)) / (1 + sin(theta_omega_max)^2);
else
   R_omega = 0;
   eta_omega = 0;
end

%% Final total acceleration
f_LT = target_flag(1) * (eta_a > tol) * (1 - (a == a_target)) * R_a * f_a + ...
       target_flag(2) * (eta_e > tol) * (1 - (e == e_target)) * R_e * f_e + ...
       target_flag(3) * (eta_incl > tol) * (1 - (incl == incl_target)) * R_incl * f_incl + ...
       target_flag(4) * (eta_RAAN > tol) * (1 - (RAAN == RAAN_target)) * R_RAAN * f_RAAN + ...
       target_flag(5) * (eta_omega > tol) * (1 - (omega == omega_target)) * R_omega * f_omega;

if sum(f_LT == 0) ~= 3
   f_LT = f_LT ./ norm(f_LT);
end
