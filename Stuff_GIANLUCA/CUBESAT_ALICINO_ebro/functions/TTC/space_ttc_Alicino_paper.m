function [M,P,info] = space_ttc_Alicino_paper(x,ep)

%% space_ttc: Communications system model
%
%   [M,P,info] = space_ttc(x,ep)
%
%% Inputs:
% * x: Design parameters
%       * x(1) = frequency, GHz
%       * x(2) = modulation,  
%                   PSK  = [0 1/8)
%                   BPSK = [1/8 2/8) 
%                   CFSK = [2/8 3/8) 
%                   BFSK = [3/8 4/8)
%                   FSK  = [4/8 5/8)
%                   DPSK = [5/8 6/8)
%                   QPSK = [6/8 7/8)
%                   NRZ  = [7/8 1]
%       * x(3) = antenna efficiency
%       * x(4) = antenna gain, dB
%       * x(5) = onboard loss, dB
%       * x(6) = other unmodelled losses, dB (polarization,implementation...)
%       * x(7) = mass of distribution network, kg
%       * x(8) = modulation index, rad
%       * x(9) = amplifier type, TWTA = [0 0.5), SSPA = [0.5 1]
%
% * ep: Environmental parameters
%       * ep(1) = Bit Error Rate
%       * ep(2) = data volume, bits
%       * ep(3) = ground station G/T, dB
%       * ep(4) = range, km
%       * ep(5) = elevation angle, deg
%       * ep(6) = pointing accuracy, deg
%       * ep(7) = ground station antenna diameter, m
%       * ep(8) = ground station access time, min
%       * ep(9) = link margin, dB
%       
%% Outputs:
% * M = total mass of the tt&c subsystem
% * P = power consumption of the tt&c subsystem
% * info:
%       * info.mass = mass inventory
%       * info.antenna = antenna design
%       * info.transponder = transponder design
%       * info.losses = losses
%
%% Author: Simone Alicino, 2013
%  Reviewed by Aarón del Río Bellisco 2015

% Design variables
f = x(1);
eta_t = x(3);
Gt = x(4);
Lt = -abs(x(5));
Lother = -abs(x(6));
mrfdn = x(7);
beta = x(8);
if x(9) < 0.5
    amplifier = 'TWTA';
else
    amplifier = 'SSPA';
end

if x(2) >= 0 && x(2) < 1/8
    modulation = 'PSK';
elseif x(2) >= 1/8 && x(2) < 2/8
    modulation = 'BPSK';
elseif x(2) >= 2/8 && x(2) < 3/8
    modulation = 'CFSK';
elseif x(2) >= 3/8 && x(2) < 4/8
    modulation = 'BFSK';
elseif x(2) >= 4/8 && x(2) < 5/8
    modulation = 'FSK';
elseif x(2) >= 5/8 && x(2) < 6/8
    modulation = 'DPSK';
elseif x(2) >= 6/8 && x(2) < 7/8
    modulation = 'QPSK';
elseif x(2) >= 7/8 && x(2) <= 1
    modulation = 'NRZ';
end 

% Environmental parameters
BER = ep(1);
V = ep(2);
GoverT0 = ep(3);
range = ep(4);
el = ep(5);
e = ep(6);
Dr = ep(7);
ta = ep(8)*60;
link_margin = ep(9);

% Hardcoded parameters
c = 299792458;
lambda = c/(1e9*f);
Re = 6378;

% Gr = GoverT0 + 10*log10(Ts0);
Gr = 10*log10(0.55*(pi*Dr/lambda)^2);
Ts0 = 10^( (Gr - GoverT0)/10 );

R = V/ta;
% R = V;

% % Slant range (actual distance between antennas)
% gamma = acosd( Re/(Re+h) );
% range = sqrt( Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cosd(gamma - el) );

% High gain antenna
[M_antenna, antenna] = antenna_design(lambda, eta_t, Gt);
Gt = antenna.gain;

% Omni-directional antenna
M_lga = 2*0.2;

% Free space loss
% Ls = 20*log10( lambda/(4*pi*range) );
Ls = -(92.44 + 20*log10(range) + 20*log10(f));

% Atmospheric attenuation
La = -atmloss(f,el);
% Increase in antenna temperature due to atmosphere
ev = [1 90];
fv = [0.1 0.5 0.8 1 2 3 4 5 6 7 8 9 10 15 20];
at = [8000 400;
    200 14;
    90 4;
    80 2.5;
    70 2;
    100 2;
    120 2;
    130 2;
    120 2;
    120 2.5;
    130 3;
    140 3.5;
    150 3.5
    220 5;
    250 10;];
Ta = interp2(ev, fv, at, el, f,'linear');

% Rain attenutation
% Lr0 = [10 0 2 4 8 40;
%        20 0 0 2 5 15;
%        45 0 0 2 3  6];
% Lr = -interp2(log10([2 8 10 20 40]),Lr0(:,1),Lr0(:,2:6),log10(f),el,'linear');
% % Lr = -rainloss(f,el);
% % Increase in antenna temperature due to rain
% Tr0 = 290; % rain temperature, K
% Tr = Tr0*(1 - 10^(Lr/10));
Lr = 0;
Tr = 0;

GoverT = Gr - 10*log10(Ts0 + Ta + Tr);

% Pointing loss
Lp = -12*(e/antenna.beamwidth)^2;

% Modulation loss (data, square wave)
Lm = 20*log10(sin(beta));

% Total losses
Ltotal = La + Lr + Lp + Lm + Lother;

% Link budget
EbNo = BER2EbNo(BER, modulation, beta);
PtdB = EbNo - Gt - Lt - Ls - Ltotal - GoverT + 10*log10(R) - 228.6 + link_margin;
Pt = 10^(PtdB/10);

% Transponder
[M_transponder, P_transponder, transponder] = transponder_design(Pt, amplifier);

% RFDU
M_rfdn = rfdn(mrfdn);

% Total mass and power
M = M_antenna + M_transponder + M_rfdn + M_lga;
P = P_transponder;

mass.antenna = M_antenna;
mass.transponder = M_transponder;
mass.rfdn = M_rfdn;
mass.total = M;
info.mass = mass;
info.antenna = antenna;
info.transponder = transponder;
losses.free_space = Ls;
losses.atmosphere = La;
losses.rain = Lr;
losses.pointing = Lp;
losses.other = Lother;
info.losses = losses;

end

function [M_antenna, antenna] = antenna_design(lambda, eta, G)

G1 = 10^(G/10);
D1 = lambda/pi*sqrt(G1/eta);
theta = 41253*eta/G1;

if G < 10
    % Patch antenna
    D = D1;
    antenna.type = 'patch';
    A = pi*D^2/4;
    rho_copper = 8940; % kg/m^3
    rho_diel = 2000;   % kg/m^3
    M_antenna = A*(0.0005*rho_copper + 0.0015*rho_diel);
elseif G < 20 || D1 < 0.4
    % Horn antenna
    antenna.type = 'horn';
    D = lambda*10^((G-7)/20);
%     theta = 72*lambda/D;
    L = D^2/3/lambda;
    H = sqrt(L^2 - D^2/4);
    S = pi*D*H/2;
    rho = 15; % surface density, kg/m^2
    M_antenna = S*rho;
    antenna.length = H;
else
    % Parabolic antenna
    D = D1;
    antenna.type = 'parabola';
    M_antenna = 2.89*D^2 + 6.11*D - 2.59;
%     rho = 10; % surface density, kg/m^2
%     M_antenna = pi*D^2*rho/4;
end

antenna.diameter = D;
antenna.gain = G;
antenna.efficiency = eta;
antenna.beamwidth = theta;
antenna.mass = M_antenna;

end

function [M_transponder, P_transponder, transponder] = transponder_design(Pt, amplifier)

% Amplified transmitter mass and power consumptions, from SMAD

transponder.power_output = Pt;

switch amplifier
    case 'SSPA'
    % Solid state
    transponder.amplifier = 'solid state';
    x1 = [1.4 2   4   6   8   10  20  30  40  50  60  70  80];
    y1 = [1.2 1.2 1.2 1.2 1.2 1.2 1.3 1.5 1.7 1.9 2.2 2.8 3.5];
    M = interp1(log10(x1),log10(y1),log10(Pt),'linear','extrap');
    transponder.mass = 10^M;
    
    x2 = [0.1 0.5  2 4   6  8 10 20 30 40 50   60  70];
    y2 = [9   12  20 28 30 32 35 50 70 90 130 180 200];
    P = interp1(log10(x2),log10(y2),log10(Pt),'linear','extrap');
    transponder.power_input = 10^P;
    
    case 'TWTA'
    % TWTA
    transponder.amplifier = 'twta';
    x1 = [0.5 0.7 0.9 2   4   6 8   10  20 30 40   50  60  70 80];
    y1 = [1.8 1.9 2   2.3 2.8 3 3.3 3.6 4  5  5.8  6.2 6.5  7  8];
    M = interp1(log10(x1),log10(y1),log10(Pt),'linear','extrap');
    transponder.mass = 10^M;
    
    x2 = [0.1 0.5  1  2   4   6   8  10  20 30 40   50   60  70  80  90];
    y2 = [8   9   10 14  16  18  21  24  32 45 60   80  105 130 155 180];
    P = interp1(log10(x2),log10(y2),log10(Pt),'linear','extrap');
    transponder.power_input = 10^P;
end

M_transponder = transponder.mass*2; % x2 redundancy
P_transponder = transponder.power_input;

%Tranponder sizing [ARB]
density_transponder=1258.58731; %obtained from different models and manufacturers
transponder.volume = M_transponder/density_transponder;
transponder.side=(transponder.volume)^(1/3); %transponder sized as a cube


end

function M_rfdn = rfdn(mrfdn)

M_rfdn = mrfdn;

end

function EbNo = BER2EbNo(BER, modulation, beta)

switch upper(modulation)
    
    case {'PSK','BPSK'}
        EbNo = (erfcinv(2*BER))^2;
        
    case 'CFSK'
        EbNo = 1/0.6*(erfcinv(2*BER))^2;
        
    case 'BFSK'
        EbNo = sqrt(2)*(erfcinv(2*BER));
        
    case 'FSK'
        EbNo = -2*log(2*BER);
        
    case 'DPSK'
        EbNo = -log(2*BER);
        
    case 'QPSK' % BER same as BPSK, but needs twice as much power. CHECK
        EsNo = 2*(erfcinv(2*BER))^2;
        EbNo = 0.5*EsNo;
        
    case 'NRZ'
        k = 0.5*( 1 - sin(2*pi*beta)/(2*pi*beta) );
        EbNo = 1/k*(erfcinv(2*BER));
        
end

end

function A = atmloss(f,el)

% From Annex 2, ITU-R P.676-9
% NOTE: 5deg <= el <= 90deg; k0 valid for f <= 54GHz

% Standard constants
p = 1013;   % total air pressure, hPa
t = 15;     % temperature, C
rho = 7.5;  % water vapour density, g/m3

% Parameters
rp = p/1013;
rt = 288/(273 + t);
phi = @(rp,rt,a,b,c,d) rp^a * rt^b * exp( c*(1-rp) + d*(1-rt) );
x1 = phi(rp,rt,0.0717,-1.8132,0.0156,-1.6515);
x2 = phi(rp,rt,0.5146,-4.6368,-0.1921,-5.7416);
x3 = phi(rp,rt,0.3414,-6.5851,0.2130,-8.5854);
n1 = 0.955*rp*rt^0.68 + 0.006*rho;
n2 = 0.735*rp*rt^0.5 + 0.0353*rt^4*rho;
g = @(f,fi) 1 + ( (f-fi)/(f+fi) )^2;
p = @(a,b,c,d,n) a*n*exp( b*(1-rt) )/( (f-c)^2 + d*n^2 );
t1 = 4.64*exp(-( (f-59.7)/(2.87+12.4*exp(-7.9*rp)) )^2)/(1 + 0.066*rp^-2.3);
t2 = 0.14*exp(2.12*rp)/( (f-118.75)^2 + 0.031*exp(2.2*rp) );
t3 = 0.0114*f/(1 + 0.14*rp^-2.6)*(-0.0247+0.0001*f+1.61e-6*f^2)/(1-0.0169*f+4.1e-5*f^2+3.2e-7*f^3);
sw = 1.013/( 1 + exp(-8.6*(rp-0.57)) );

% Specific attenuations due to dry air and water vapour, dB/km
k0 = ( (7.2*rt^2.8)/(f^2 + 0.34*rp^2*rt^1.6) + (0.62*x3)/((54-f)^(1.16*x1) + 0.83*x2) )*f^2*rp^2*1e-3;
kw = ( p(3.98,2.23,22.235,9.42,n1)*g(f,22) + p(11.96,0.7,183.31,11.14,n1)...
    + p(0.081,6.44,321.226,6.29,n1) + p(3.66,1.6,325.153,9.22,n1) + p(25.37,1.09,380,0,n1)...
    + p(17.4,1.46,448,0,n1) + p(844.6,0.17,557,0,n1)*g(f,557) + p(290,0.41,752,0,n1)*g(f,752)...
    + p(8.3328e4,0.99,1780,0,n2)*g(f,1780) )*f^2*rt^2.5*rho*1e-4;

% Equivalent height for dry air and water vapour, km
h0 = 6.1*(1 + t1 + t2 + t3)/(1 + 0.17*rp^-1.1);
hw = 1.66*(1 + 1.39*sw/((f-22.235)^2 + 2.56*sw) + 3.37*sw/((f-183.31)^2 + 4.69*sw) + 1.58*sw/((f-325.1)^2 + 2.89*sw) );

% Path attenuation, dB
A0 = k0*h0;
Aw = kw*hw;
A = (A0 + Aw)/sind(el);

end
