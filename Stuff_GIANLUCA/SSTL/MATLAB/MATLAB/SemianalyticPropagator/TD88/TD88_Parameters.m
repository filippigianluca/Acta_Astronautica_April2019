function [fo,fx,ko,g1_4,gB,delta,Tau,cj,zj,Bj]=TD88_Parameters(X,solardata,constant,d,t,c1,p,H,S,m)

% orbital elements
a=X(1);             % semimajor axis
e=X(2);             % eccentricity
I=X(3);             % inclination
w=X(4);             % argurment of perigee
CapitalOmega=X(5);  % RAAN
f=X(6);             % true anomaly

%%%% Parameters %%%%%

% Parameter accountable for variable solar and gemagnetic activity

% %%%%%% Change this later
% fm=(100-60)/160;
% 
% fo=c1(2)+fm;
% 
% fx=1+c1(1)*(100-100);
% 
% ko=1+c1(3)*(6-3);



fm=(solardata.Fb(t)-60)/160;

fo=c1(2)+fm;

fx=1+c1(1)*(solardata.Fx(t)-solardata.Fb(t));

ko=1+c1(3)*(solardata.kp(t)-3);





% g(1-4) functions
g1_4=[1, ((fm/2)+c1(4)), (1+c1(5)*fm)*sin(2*(d-p(4))), (1+c1(6)*fm)*sin(d-p(5))];


% gBar(5-7) functions
gB=[0,0,0,0,sin(d-p(3)), (1+c1(7)*fm), (1+c1(8)*fm)];


% the inverse of the ballistics coefficent, delta
delta=constant.Cd*S/m;


% tau
Tau=constant.omegat * sqrt(1 - (e^2)) * cos(I) / (sqrt(constant.mu / (a^3)));


% cj
cj=[1, 1/2, 1/3]'.*(constant.eta * constant.Rearth * (sin(I)^2) / (2 * H));


% zj
zj=[1, 1/2, 1/3]'.*(a * e / H);


% Bj
Bj=exp(((constant.hB + constant.Rearth - a)./([1, 2, 3]'*H)) - cj);
end