addpath(genpath('../../spaceart_toolbox'))

% clear
close all
clc

tic

one_day = 1 / 365.25;

% -------------------------------------------------------------------------
% Initial date of the propagation
% -------------------------------------------------------------------------
StartDate.day=19;
StartDate.month=4;
StartDate.year=1995;
to=StartDate.year+((StartDate.month-1)*28+StartDate.day)*one_day; % 0.00274 repesents on day in a year

%%% Marilena
to = StartDate.year +( date2mjd2000([StartDate.year, StartDate.month, StartDate.day, 0, 0, 0]) - ...
                      date2mjd2000([StartDate.year, 1, 1, 0, 0, 0]) +1) * one_day;

%%%%

% to=datetime('17-05-1975','InputFormat','dd-MM-yyyy');
% to=juliandate(Idate)-juliandate(datetime('01-10-1957 23:59:59','InputFormat','dd-MM-yyyy HH:mm:ss'));

% End Date
End.day=19;
End.month=5;
End.year=1995;
End.Date=End.year+((End.month-1)*28+End.day)*one_day;
% -------------------------------------------------------------------------
% Initial osculating keplerian elements (semimajor axis, eccentricity,
% inclination, right ascension, perigee argument, true anomaly)
% -------------------------------------------------------------------------
a=6768137; % meters
e=0.0016;
I= 51.645 * pi/180;
CapitalOmega=0 * pi/180;
w=90 * pi/180;
f=0 * pi/180;
X = [a e I w CapitalOmega f];
% -------------------------------------------------------------------------
% Mass of the spacecraft [kg]
% -------------------------------------------------------------------------
m = 20.63;

% -------------------------------------------------------------------------
% Cross section of the spacecraft [km^2]
% -------------------------------------------------------------------------
radius=(0.215/2); % meters 
S=(radius^2)*pi;

% -------------------------------------------------------------------------
% final altitude of propagation (km)
% -------------------------------------------------------------------------
% end_altitude=200; % km
end_altitude=150000; % meters

%% ------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------

% drag coefficent (vallado page 551)
constant.Cd=2.2;


% the speed of the rotation of the earth (rad/s) (page132 vallado)
constant.omegat=7.292115e-5; 


% Radius of earth (m) (vallado last page)
constant.Rearth=6378137; % meters

% mass of earth (kg) (vallado last page)
constant.Mearth=5.9733328e24;

% Gravitational constant (m^3*kg^-1*s^-2)
constant.G=6.673e-11; % meters

% gravitational parameter (m^3*s^-2) (vallado page 130)
constant.mu=constant.G*constant.Mearth;

% constant eta
constant.eta=1/298;

% hBar (m)
constant.hB=120000; % meters

% J2 and J3
constant.J2=0.0010826267;
constant.J3 = -2.5327e-6;

constant.one_day = one_day;


%% ------------------------------------------------------------------------
% Call Scripts
% -------------------------------------------------------------------------
% Call solar data
TD88_SolarData; % Input data of solar activity

%% ------------------------------------------------------------------------
% While loop
% -------------------------------------------------------------------------


t=to;
Results.a=[];
Results.e=[];
Results.I=[];
Results.w=[];

B = -0.5 * (constant.Rearth/a)  * (constant.J3/constant.J2) * sin(I);
h0 = e * sin(w);
k0 = e * cos(w);
aph = atan2(h0 - B,k0);
A = sqrt(k0^2 + (h0-B)^2);

h = e * sin(w);
k = e * cos(w);

% while ((X(1)-constant.Rearth)>end_altitude) && (t<End.Date)
while ((X(1)-constant.Rearth)>end_altitude) && (t<to + 30/365)
    
Results.a=[Results.a; a];
Results.e=[Results.e; e];
Results.I=[Results.I; I];
Results.w=[Results.w; w];
t=t+one_day; % one day time step
[Da,De,DI,eDw,DP1,DP2,Dw, Draan]=ChangeInOrbitalElements(t,X,m,S,solardata,constant);

%%%% Marilena 2
% X(1) = a + Da;
% a = X(1);
% X(2) = e + De;
% X(3) = X(3) + DI;
% X(4) = X(4) + Dw;
%%%%end Marilena 2


[e,w] = Geop_effect_low_e(X, constant, aph, A, t, to, one_day);

h = h + DP1;
k = k + DP2;

% Marilena 2 - comment
a = a + Da;
% end Marilena 2

% e = sqrt(h^2 + k^2);
e = e + De;


% %%% Marilena 1
w = w + Dw;
% %%% end Marilena 1


% Marilena 2 - comment
I = I + DI;
% end Marilena 2


% w = w + eDw;
% w = atan2(h , k);
% CapitalOmega = CapitalOmega + Draan;
X = [a e I w CapitalOmega f];

% hold on;plot(t,mod(w, 2*pi)*180/pi,'.','color','b');pause(0.005)
end

%% Plots

figure
subplot(2,2,1)
plot(to:one_day:t,Results.a./1000,'color','b')
title('Semimajor axis, a')
xlabel('Time (years)')
ylabel('Semimajor axis, a (km)')
% axis([1995 2000 6500 6770])
grid on
set(gca,'fontsize',20,'fontweight','demi')

subplot(2,2,2)
plot(to:one_day:t,Results.e,'color','b')
title('Eccentricity, e')
xlabel('Time (years)')
ylabel('Eccentricity, e')
% axis([1995 2000 0.00005 0.0017])
set(gca,'fontsize',20,'fontweight','demi')
grid on

subplot(2,2,3)
plot(to:one_day:t,Results.I.*(180/pi),'color','b')
title('Inclination, I')
xlabel('Time (years)')
ylabel('Inclination, I (degrees)')
% axis([1995 2000 51.6 51.65])
set(gca,'fontsize',20,'fontweight','demi')
grid on

subplot(2,2,4)
plot(to:one_day:t,Results.w.*(180/pi),'color','b')
title('Argument of Perigee, w')
xlabel('Time (years)')
ylabel('Argument of Perigee, w (degrees)')
% axis([1995 2000 0 180])
set(gca,'fontsize',20,'fontweight','demi')
grid on

decay_days=round((t-to)/0.00274);
percentage_error=(1525 - decay_days) * 100 / 1525;

disp('Computation complete')
disp('Number of days for decay to complete =')
disp(decay_days)
disp('Percentage error from actual value =')
disp(percentage_error)