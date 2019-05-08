function [DeltaV, a, e, omega, m, acc, thrust] = Fernandes_a_e_omega(initial_orbit, final_orbit, t, engine, m0, constants)

% Input: initial_orbit -> structure containing a, e and omega of initial
%                         orbit
%        final_orbit -> structure containing a, e and omega of the final
%                       orbit

% References: "Optimal Low-Thrust Transfers between coplanar orbits with
% small eccentricities", Sandro da Silva Fernandes at al
% "A first order analytical theory for optimal low thrust limited power
% transfers between arbitrary elliptical coplanar orbits", Sandro da Silva
% Fernandes et al

% Marilena Di Carlo, 2015

% Final time
T = t(end);

% Initial orbit
a0 = initial_orbit.a;
e0 = initial_orbit.e;
omega0 = initial_orbit.omega;

% Final orbit
af = final_orbit.a;
ef = final_orbit.e;
omegaf = final_orbit.omega;

% Initial h and k
h0 = e0 * cos(omega0);
k0 = e0 * sin(omega0);

% Final h and k
hf = ef * cos(omegaf);
kf = ef * sin(omegaf);

% Final variations of h and k
Delta_hf = hf - h0;
Delta_kf = kf - k0;

% alpha f
alpha_f = af / a0;

% psi f (Equation 42 Ref. 1)
psi_f = sqrt( (Delta_hf^2 + Delta_kf^2) / 5);

% K0 (not to be confused with k = e sin(omega)
K0 = atan2( sin(sqrt(2) * psi_f) , (sqrt(alpha_f) - cos(sqrt(2) * psi_f)));

% Equation 44 Ref. 1 applied at final time
uf_v0 = sqrt(1 - 2 * sqrt(a0/af) * cos( sqrt(2) * psi_f) + (a0/af));

v0 = sqrt(constants.mu/a0);

% E
E = v0^2 / (2 * T^2) * (uf_v0)^2;

% C
C = sqrt( (4 * constants.mu * E * (sin(K0))^2) / (5 * a0));


%% Initial conditions adjoint variables

% Initial condition adjoint variable a
pa0 = sqrt(E/2) * (v0/a0) * cos(K0);

% Initial condition adjoint variable h
ph0 = Delta_hf * C / (psi_f * sqrt(5));

% Initial condition adjoint variable k
pk0 = Delta_kf * C / (psi_f * sqrt(5));


%% Semimajor axis

% Equation 29 ref.1
a = a0 ./ (1 + 4 * a0 / constants.mu * (0.5 * E * t.^2 - a0 * pa0 * t));

% Equation 32 ref 1
pa = sqrt( (a0./a).^3 * pa0^2 + 5/8 * C^2 * (a0./a.^3 - 1./a.^2));


%% h

% Equation 33 ref.1
ph = ph0;   

% Equation 30 ref.1
h = h0 + sqrt(2.5) * ph / C * ( atan2(sqrt(4 * constants.mu * E - 5 * C^2 * a0), sqrt(5 * C^2 * a0)) - ...
                                atan2(sqrt(4 * constants.mu * E - 5 * C^2 * a) , sqrt(5 * C^2 * a)) );
                         
                            
%% k

% Equation 34 ref. 1
pk = pk0;

% Equation 31 ref. 1
k = k0 + sqrt(5/2) * pk / C * ( atan2( sqrt(4 * constants.mu * E - 5 * C^2 * a0), sqrt(5 * C^2 * a0)) - ...
                                atan2( sqrt(4 * constants.mu * E - 5 * C^2 * a), sqrt(5 * C^2 * a)) );
                            
                        
                            
%% Eccentricity

e = sqrt(h.^2 + k.^2);

%% Argument of the perigee

omega = atan2(k ./e, h./e);


%% J [km^3/s^3]

Delta_h = h - h0;
Delta_k = k - k0;

% Equation 42 ref 1
psi = sqrt(( Delta_h.^2 + Delta_k.^2 ) / 5);

% Equation 45 ref. 1
J  = (v0^2 ./ ( 2 * t)) .* (1 - 2 * sqrt(a0./a) .* cos( sqrt(2) * psi) + (a0 ./ a));


%% Mass [kg]

% The following equation is obtianed combining equation 45 and 1
m = m0 * engine.P_max ./ (J * m0 + engine.P_max);


%% Acceleration [km/s^2]

acc = engine.P_max ./ (constants.g0 * engine.Isp * m);
 

%% Thrust - PROBABLY THIS IS WRONG!!!!!

% Thrust in kg km/s^2
thrust = m .*acc;

% Thrust in N
thrust = thrust * 1000;
 
%% DeltaV [km/s]

% DeltaV =integral(@(tt)deltaV_Fernandes(tt, a0, k0, h0, pa0, pk0, ph0, E, C, mu,P_max,Isp,g0,m0), 0, t(end));
 
DeltaV = engine.Isp * constants.g0 * log(m0 / m(end));
 