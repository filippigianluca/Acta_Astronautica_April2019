function rho = TD88(x, Fx, Fb, Kp, d, alpha_Sun)

% Constants - write outside
R = 6378.136e3;
epsilon = 1/298;


% Constants of TF88 model - write outside
h_bar = 120e3;
H = 29e3;
a = [0.007, 0.2875, 0.04762, 0.0471, 7, 7, 0.3333, 15];
p = [0, 0, -29.41, -263, 263, 8.0913, 10.0813];
Knj = [ .766373e-8,   .16573e-9,   .3871e-10; ...
       -.440149e-8,   .33428e-9,   .9352e-10; ...
        .118107e-9,  -.14781e-9,  -.1518e-11; ...
       -.159664e-10, -.64670e-11, -.2050e-11; ...
       -.240755e-9,  -.13985e-10, -.3059e-11; ...
        .643785e-10,  .13618e-9,   .3517e-10; ...
        .744666e-11,  .45416e-11,  .2080e-11];
Kn = [0.296815e-14; 2.81456e-13; -.123300e-13; -.114892e-16; -.390065e-15; ...
       .742439e-14; -.341594e-15];

   
% Equinoctial and classical elements

% Equinoctial elements
sma  = x(1);
P1 = x(2);
P2 = x(3);
Q1 = x(4);
Q2 = x(5);
l  = x(6);

% Classical elements
e          = sqrt(P1^2 + P2^2);
incl       = 2 * atan( sqrt(Q1^2 + Q2^2) );
RAAN       = atan2(Q1, Q2);
omega_RAAN = atan2(P1, P2);
omega      = mod(omega_RAAN - RAAN, 2*pi);

% Solve Kepler equation in the eccentirc longitude K = Omega + omega + E
K = kepler_K(e, l, P1, P2);
E = mod(K - omega_RAAN, 2*pi);

zeta = mod(omega + E, 2*pi);
sin_phi = sin(incl) * sin(zeta);
F_ = 1- e * cos(E);
h = sma * F_ - R * (1 - epsilon * sin_phi^2);



% Useful quantities
cos_t__cos_phi = - cos(incl) * sin(zeta) * sin(alpha_Sun - RAAN) + ...
              - cos(zeta) * cos(alpha_Sun - RAAN);
          
sin_t__cos_phi = - cos(incl) * sin(zeta) * cos(alpha_Sun - RAAN) + ...
                cos(zeta) * sin(alpha_Sun - RAAN);

fx = 1 + a(1) * (Fx - Fb);
fm = (Fb - 60) / 160;
f0 = a(2) + fm;
k0 = 1 + a(3) * (Kp - 3);

hn_exp = exp( (h_bar - h) / H);
hn_exp1 = hn_exp * Knj;
hn = Kn + sum(hn_exp1,2);

g(1,1) = 1;
g(2,1) = fm/2 + a(4);
g(3,1) = (1 + a(6) * fm) * sin(2*(d-p(3)));
g(4,1) = (1 + a(5) * fm) * sin(d-p(4));
g(5,1) = sin(d-p(5)) * sin(incl) * sin(omega+E);
g(6,1) = (1 + a(7) * fm) * (cos(p(6)) * sin_t__cos_phi - sin(p(6)) * cos_t__cos_phi);
g(7,1) = (1 + a(8) * fm) * (2 * cos(2 * p(7)) * sin_t__cos_phi  * cos_t__cos_phi + ...
                            - sin(2 * p(7)) * (sin_t__cos_phi)^2);

rho = fx * f0 * k0 * dot(g, hn);

% keyboard