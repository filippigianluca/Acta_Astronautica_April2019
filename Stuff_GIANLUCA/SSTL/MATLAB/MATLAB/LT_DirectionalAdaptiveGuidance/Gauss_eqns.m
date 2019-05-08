function dOE_dt = Gauss_eqns(t, x, DAG_parameters, engine, constants)


T = engine.T;
Isp = engine.Isp;



% Orbital elements
a    = x(1);
e    = x(2);
incl = x(3);
RAAN = x(4);
omega = x(5);
theta = x(6);

% Mass
m     = x(7);

% Acceleration (magnitude)
f = T / m;

% Acceleration vector
f_LT = DAG_LT_acc(x, DAG_parameters, constants);
f_LT = f * f_LT;

fR = f_LT(1);
fC = f_LT(2);
fN = f_LT(3);

p = a * (1 - e^2);
h = sqrt(constants.mu * p);
r = p / (1 + e * cos(theta));

% Gauss equations
da_dt = 2*a^2 / h * (e * sin(theta) * fR + p/r * fC);

de_dt = 1/h * (p * sin(theta) * fR + ...
        ( (p+r) * cos(theta) + r * e) * fC);
    
dincl_dt = r/h * cos(theta) * fN;

dRAAN_dt = r/h * sin(theta) / sin(incl) * fN;

domega_dt = 1 / (h*e) * (- p * cos(theta) * fR + ...
            (p+r) * sin(theta) * fC) + ...
           - r/h * sin(theta) * cot(incl) * fN;
       
dtheta_dt = h/r^2 + 1/(e*h) * (p * cos(theta) * fR - ...
            (p+r) * sin(theta) * fC);
      
% Propellant mass reduces only when engine is on
dm_dt = - (sum(f_LT==0)~=3) * T / (Isp * constants.g0);        


dOE_dt = [da_dt; de_dt; dincl_dt; dRAAN_dt; domega_dt; dtheta_dt; dm_dt];