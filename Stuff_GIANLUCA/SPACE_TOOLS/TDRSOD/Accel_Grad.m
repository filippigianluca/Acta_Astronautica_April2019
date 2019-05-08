%--------------------------------------------------------------------------
%
% Accel_Grad: Computes the acceleration and gradient of an Earth orbiting
%             satellite due to the Earth's gravity field (low degree and
% 			  order)
%
% Inputs:
%   Mjd_UTC     Modified Julian Date (UTC)
%   r           Satellite position vector in the ICRF/EME2000 system
%   v           Satellite velocity vector in the ICRF/EME2000 system
%   Area        Cross-section 
%   mass        Spacecraft mass
% 
% Outputs:
%   a           Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%   G           Gradient (G=da/dr) in the ICRF/EME2000 system
%   dadCD       Partials of acceleration w.r.t. to drag coefficient
%   dadCR       Partials of acceleration w.r.t. to solrad coefficient
%
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function [a, G, dadCD, dadCR] = Accel_Grad(x, Y)

global const AuxParam eopdata

MJD_UTC = AuxParam.Mjd_UTC+x/86400;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_TT = MJD_UTC+TT_UTC/86400;
MJD_UT1 = MJD_UTC+UT1_UTC/86400;

P = PrecMatrix(const.MJD_J2000,MJD_TT);
N = NutMatrix(MJD_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(MJD_UT1) * T;

% Difference between ephemeris time and universal time
JD     = MJD_UTC+2400000.5;
[year, month, day, hour, minute, sec] = invjday(JD);
temp   = JD-2415019.5;
leapyrs= floor((year-1901)*0.25);
days   = temp-((year-1900)*365+leapyrs);
ET_UT  = ETminUT(year+days/365.25);
MJD_ET = MJD_UTC+ET_UT/86400;
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_ET);

% MJD_TDB = Mjday_TDB(TT);
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_TDB);

% Acceleration due to the Earth's gravity field
if (AuxParam.SolidEarthTides || AuxParam.OceanTides)
	a = AccelHarmonic_ElasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole,y_pole);
    % a = AccelHarmonic_AnelasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole,y_pole);
else
    a = AccelHarmonic(Y(1:3),E,AuxParam.n,AuxParam.m);
end

% Gradient
G = EarthHarmGrad(Y(1:3),E,AuxParam.n,AuxParam.m);

% Atmospheric density
dens = nrlmsise00(MJD_UTC,T*Y(1:3),UT1_UTC,TT_UTC);

% Drag coefficient partials
Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-const.MJD_J2000)/36525 ); % [rad/s]
% Omega = const.omega_Earth-0.843994809*1e-9*LOD; % IERS [rad/s]
dadCD = AccelDrag(dens,Y(1:3),Y(4:6),T,AuxParam.area_drag,AuxParam.mass,1,Omega);

% Radiation pressure coefficient partials
dadCR = AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB,...
        AuxParam.area_solar,AuxParam.mass,1,const.P_Sol,const.AU,'geometrical');

