%--------------------------------------------------------------------------
%
% Range_4W.m: Computes the TDRS 4-way range
%
% Inputs:
%   Mjd_UTC     Ground-received time of measurement
%   Trj_User    User satellite trajectory integration object
%   Trj_TDRS    Relay satellite trajectory integration object
%   Sta         Ground station
%
% Outputs:
%   rho         4-way range
%   drho_dUser  Partials w.r.t. user satellite epoch state and parameters
%   drho_dTDRS  Partials w.r.t. TDRS satellite epoch state and parameters
%
% Last modified:   2018/01/04   M. Mahooti
%
%--------------------------------------------------------------------------
function [rho,drho_dUser,drho_dTDRS] = Range_4W(Mjd_UTC,Trj_User,Trj_TDRS,Sta)

global const eopdata

i_max = 2;                  % Order of light-time correction
n_var = 8;

PhiS = zeros(6,n_var);      % transition and sensitivity matrix

% Ground received time (Terrestrial Time)
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;
Mjd_TT = Mjd_UTC + TT_UTC/86400;

% Precession, nutation
PN = (NutMatrix(Mjd_TT) * PrecMatrix(const.MJD_J2000,Mjd_TT))';

% Station (pseudo Earth-fixed)
R_Sta = PoleMatrix(x_pole,y_pole)' * Position(Sta.lon,Sta.lat,Sta.h)';

% Downleg TRDS -> station
tau1 = 0;
r_Sta = PN * GHAMatrix(Mjd_UT1)' * R_Sta;

for i=1:i_max+1
    t = ( -tau1 )*86400;
    y_TDRS = DEInteg(@Accel,0,t,1e-13,1e-6,6,Trj_TDRS.y);
    r_TDRS = y_TDRS(1:3);
    r1     = r_TDRS-r_Sta;
    rho1   = norm(r1);
    tau1   = (rho1/const.c_light)/86400;
end

% Downleg User -> TDRS
tau2 = 0;
for i=1:i_max+1
    t = ( -tau1-tau2 )*86400;
    y_User = DEInteg(@Accel,0,t,1e-13,1e-6,6,Trj_User.y);
    r_User = y_User(1:3);
    r2     = r_User-r_TDRS;
    rho2   = norm(r2);
    tau2   = (rho2/const.c_light)/86400;
end

% Upleg TDRS -> User
tau3 = 0;
for i=1:i_max+1
    t = ( -tau1-tau2-tau3 )*86400;
    y_TDRS = DEInteg(@Accel,0,t,1e-13,1e-6,6,Trj_TDRS.y);
    r_TDRS = y_TDRS(1:3);
    r3     = r_TDRS-r_User;
    rho3   = norm(r3);
    tau3   = (rho3/const.c_light)/86400;
end

% Upleg station -> TDRS
tau4 = 0;
for i=1:i_max+1
    r_Sta  = PN * (GHAMatrix(Mjd_UT1-tau1-tau2-tau3-tau4))' * R_Sta;
    r4     = r_Sta-r_TDRS;
    rho4   = norm(r4);
    tau4   = (rho3/const.c_light)/86400;
end

% Range
rho = 0.5 * ( rho1 + rho2 + rho3 + rho4 );

% Range partials w.r.t. parameters of user satellite 
t = ( -tau1-tau2 )*86400;
Y_User = DEInteg(@VarEqn_,0,t,1e-13,1e-6,7*n_var,Trj_User.Y);
for j=1:n_var
    PhiS(:,j) = Y_User(6*j+1:6*j+6);
end
drho_dUser =  0.5 * ( r2/rho2 - r3/rho3 )' * PhiS(1:3,1:8);

t = ( -tau1 )*86400;
Y_TDRS = DEInteg(@VarEqn_,0,t,1e-13,1e-6,7*n_var,Trj_TDRS.Y);
for j=1:n_var
    PhiS(:,j) = Y_TDRS(6*j+1:6*j+6);
end
drho_dTDRS =  0.5 * ( r1/rho1 - r2/rho2 )' * PhiS(1:3,1:8);

t = ( -tau1-tau2-tau3 )*86400;
Y_TDRS = DEInteg(@VarEqn_,0,t,1e-13,1e-6,7*n_var,Trj_TDRS.Y);
for j=1:n_var
    PhiS(:,j) = Y_TDRS(6*j+1:6*j+6);
end

drho_dTDRS = drho_dTDRS + 0.5 * ( r3/rho3 - r4/rho4 )' * PhiS(1:3,1:8);

