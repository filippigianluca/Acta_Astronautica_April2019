global  AuxParam
global PC
load DE430Coeff.mat
PC = DE430Coeff;

AuxParam = struct ('Mjd_UTC',0,'area_drag',0,'area_solar',0,'mass',0,'Cr',0,...
                   'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                   'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

 MJD_UTC = AuxParam.Mjd_UTC+t/86400;              
               
               
% % Initialization
% t = (Mjd_UTC-AuxParam.Mjd_UTC)*86400;
                
               
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


%        Moon wrt Earth          pccor      rpc
%        Earth wrt Sun           ccor       rc
%        Moon wrt Sun            pscor      rps   
%        Satellite wrt Earth     sbcor      rsb  
%        Satellite wrt Sun       bcor       rb 
%        Satellite wrt Moon      sbpcor     rsbp
pccor = r_Moon;
ccor = r_Earth-r_SunSSB;
pscor = r_Moon-r_Sun;
sbcor = r;
bcor = r-r_Sun;
sbpcor = r-r_Moon;


[nu,~] = Shadow(pccor,ccor,pscor,sbcor,bcor,sbpcor);