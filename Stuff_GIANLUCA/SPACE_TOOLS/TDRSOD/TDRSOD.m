%--------------------------------------------------------------------------
%
% TDRSOD: Orbit Determination from Tracking and Data Relay Satellite
%         measurements
%
% Last modified:  2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
clc
clear
format long g
profile on

global const CS Cnm Snm AuxParam eopdata swdata R d PC

SAT_Const
EGM96
load DE430Coeff.mat
PC = DE430Coeff;

Cnm = zeros(71,71);
Snm = zeros(71,71);
fid = fopen('egm96','r');
for n=0:70
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end

AuxParam = struct ('Mjd_UTC',0,'area_drag',0,'area_solar',0,'mass',0,'Cr',0,...
                   'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                   'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

Sta = struct ('StaNo',0,'lat',0,'lon',0,'h',0);

TrjData = struct ('y',0,'Y',0,'PhiS',0);

% TDRS tracking data record
Obs = struct ('Mjd_UTC',0,'StaNo',0,'TdrsNo',0,'Range_4w',0);

sigma_range = 10;        % TDRS measurement accuracy [m]

% Number of decimals for printout of estimation parameters
Digits = [2,2,2,5,5,5,4,4];

% Estimation parameter names
Label = ['x','y','z','D','R'];

% Read setup parameters from file
fid = fopen('TDRSOD.txt','r');

% TDRS satellite numbers
tline = fgetl(fid);
n_sat = str2num(tline(26:38));  % Number of TDRS satellites
n_sat = n_sat+1;                % Total number of satellites

TdrsNo(2) = str2num(tline(39:51));
TdrsNo(3) = str2num(tline(52:64));

TdrsNo(1) = 0;

% Station numbers
tline = fgetl(fid);
n_Sta = str2num(tline(26:38));

Sta(1).StaNo = str2num(tline(39:51));
Sta(2).StaNo = str2num(tline(52:64));

% Number of iterations
tline = fgetl(fid);
n_iterat = str2num(tline(26:38));

% Epoch
tline = fgetl(fid);
Y = str2num(tline(26:32));
M = str2num(tline(34:35));
D = str2num(tline(37:38));
h = str2num(tline(40:41));
m = str2num(tline(43:44));
s = str2num(tline(46:51));

% read Earth orientation parameters
fid1 = fopen('eop19990101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid1,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid1);

% read space weather data
fid1 = fopen('sw19990101.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid1,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid1);

Mjd0_UTC = Mjday(Y,M,D,h,m,s); % UTC
AuxParam(1).Mjd_UTC = Mjd0_UTC;
AuxParam(2).Mjd_UTC = Mjd0_UTC;
AuxParam(3).Mjd_UTC = Mjd0_UTC;

% Instantiation of vector objects
y0 = zeros(6,n_sat);                        % Epoch state (xyz,vxyz)
Sigy0 = zeros(6,n_sat);                     % Standard deviation
p = zeros(2,n_sat);                         % Parameters (Cd,Cr)
Sigp = zeros(2,n_sat);                      % Standard deviation

% Epoch state vector, s/c parameters and a priori standard deviations
% for user s/c and TDRS satellites
for i_sat=1:n_sat
    % Epoch state and standard deviation
    for i=1:6
        tline = fgetl(fid);
        y0(i,i_sat) = str2num(tline(26:38));
        Sigy0(i,i_sat) = str2num(tline(39:51));
    end
    
    % Mass, area, Cd, Cr, sig(Cd), sig(Cr)
    tline = fgetl(fid);
    AuxParam(i_sat).mass = str2num(tline(26:38));
    AuxParam(i_sat).area_drag = str2num(tline(39:51));
    AuxParam(i_sat).area_solar = AuxParam(i_sat).area_drag;
    
    tline = fgetl(fid);
    AuxParam(i_sat).Cd = str2num(tline(26:38));
    Sigp(1,i_sat) = str2num(tline(39:51));
    
    tline = fgetl(fid);
    AuxParam(i_sat).Cr = str2num(tline(26:38));
    Sigp(2,i_sat) = str2num(tline(39:51));
    
    p(1,i_sat) = AuxParam(i_sat).Cd;
    p(2,i_sat) = AuxParam(i_sat).Cr;
end

AuxParam(1).n = 20;
AuxParam(1).m = 20;
AuxParam(1).sun = 1;
AuxParam(1).moon = 1;
AuxParam(1).sRad = 1;
AuxParam(1).drag = 1;
AuxParam(1).planets = 1;
AuxParam(1).SolidEarthTides = 1;

AuxParam(2).n = 20;
AuxParam(2).m = 20;
AuxParam(2).sun = 1;
AuxParam(2).moon = 1;
AuxParam(2).sRad = 1;
AuxParam(2).drag = 1;
AuxParam(2).planets = 1;
AuxParam(2).SolidEarthTides = 1;

AuxParam(3).n = 20;
AuxParam(3).m = 20;
AuxParam(3).sun = 1;
AuxParam(3).moon = 1;
AuxParam(3).sRad = 1;
AuxParam(3).drag = 1;
AuxParam(3).planets = 1;
AuxParam(3).SolidEarthTides = 1;

Aux = AuxParam;

% Station location
Rs = zeros(3,1);
for i_Sta=1:n_Sta
    tline = fgetl(fid);
    Rs(1) = str2num(tline(26:38));
    Rs(2) = str2num(tline(39:51));
    Rs(3) = str2num(tline(52:64));
    [Sta(i_Sta).lon,Sta(i_Sta).lat,Sta(i_Sta).h] = Geodetic(Rs);
end

fclose(fid);

% Read all observations from input file to observations vector
fid = fopen('TDRSOD.dat','r');

tline = fgetl(fid);

for i=1:341
    tline = fgetl(fid);
    Y = str2num(tline(1:4));
    M = str2num(tline(6:7));
    D = str2num(tline(9:10));
    h = str2num(tline(13:14));
    m = str2num(tline(16:17));
    s = str2num(tline(19:24));
    Obs(i).Mjd_UTC = Mjday(Y,M,D,h,m,s);     % UTC
    Obs(i).StaNo = str2num(tline(28:30));
    Obs(i).TdrsNo = str2num(tline(34));
    Obs(i).Range_4w = str2num(tline(37:46));
    Obs(i).Range_4w = Obs(i).Range_4w*1000;  % [m]
end

fclose(fid);

% Allocation of dynamic variables and objects
n_est = 8*n_sat;                       % Number of estimation parameters

X_apr = zeros(n_est,1);                % Estimation parameter vector
dzdX = zeros(n_est,1);                 % Partials of modelled observations
dX = zeros(n_est,1);                   % Correction and standard deviation
SigX = zeros(n_est,1);

% A priori estimation parameter vector and covariance
for i_sat=0:n_sat-1
    offset = 8*i_sat;
    for i=1:6
        X_apr(offset+i  ) = y0(i,i_sat+1);
    end
    for i=1:2
        X_apr(offset+6+i) = p(i,i_sat+1);
    end
    for i=1:6
        SigX(offset+i  )  = Sigy0(i,i_sat+1);
    end
    for i=1:2
        SigX(offset+6+i)  = Sigp(i,i_sat+1);
    end
end

P_apr = diag(SigX)*diag(SigX); % A priori covariance

X = X_apr;

% Header
fprintf('TDRS Orbit Determination\n\n');
[r,c] = size(Obs);
fprintf('%d', c);
fprintf(' measurements read from tracking data file\n\n');

% Iteration loop
for iterat=1:n_iterat
    % Iteration header
    fprintf('Iteration ');
    fprintf('%d', iterat);
    fprintf('\n\n');
    fprintf('    Date         UTC        Sta   TDRS');
    fprintf('    obs [km]   comp [km]     o-c [m]\n');
    
    R = zeros(n_est);
    d = zeros(n_est,1);
    
    % Initialize measurement statistics
    LSQInit(n_est,X_apr-X,P_apr);
    n_obs      = 1;
    Sum_res    = 0;
    Sum_ressqr = 0;
    TrjData(1).y = y0(:,1);
    TrjData(2).y = y0(:,2);
    TrjData(3).y = y0(:,3);
    
	% State and partials
    n_var = 8;
    TrjData(1).Y = zeros(7*n_var,1);
    TrjData(2).Y = zeros(7*n_var,1);
    TrjData(3).Y = zeros(7*n_var,1);
    
    % Create combined vector from epoch state, epoch transition matrix (=1)
    % and epoch sensitivity matrix (=0)
    for i=1:6
        TrjData(1).Y(i) = TrjData(1).y(i);
        TrjData(2).Y(i) = TrjData(2).y(i);
        TrjData(3).Y(i) = TrjData(3).y(i);
        for j=1:n_var
            if (i==j)
                TrjData(1).Y(6*j+i) = 1;
                TrjData(2).Y(6*j+i) = 1;
                TrjData(3).Y(6*j+i) = 1;
            else
                TrjData(1).Y(6*j+i) = 0;
                TrjData(2).Y(6*j+i) = 0;
                TrjData(3).Y(6*j+i) = 0;
            end
        end
    end
	
    % Process measurements
    for i=1:c        
        % Integrate user and relay satellite trajectories to
        % ground-received time
        % of current measurement
        for i_sat=1:n_sat
            AuxParam = Aux(i_sat);
            [TrjData(i_sat).y,TrjData(i_sat).Y] = Trj(8,Obs(i).Mjd_UTC,TrjData(i_sat).y,TrjData(i_sat).Y);
        end
        Aux(1).Mjd_UTC = Obs(i).Mjd_UTC;
        Aux(2).Mjd_UTC = Obs(i).Mjd_UTC;
        Aux(3).Mjd_UTC = Obs(i).Mjd_UTC;
        
        % Select TDRS satellite and ground station
        I_TDRS = -1;
        for i_sat=2:n_sat 
            if (TdrsNo(i_sat)==Obs(i).TdrsNo)
                I_TDRS = i_sat-1;
                break
            end
        end
        
        if (I_TDRS<0)
            continue     % Unknown TDRS satellite; skip observation
        end
        
        % Select ground station
        I_Sta = -1;
        for i_Sta=1:n_Sta
            if (Sta(i_Sta).StaNo==Obs(i).StaNo)
                I_Sta = i_Sta-1;
                break
            end
        end
        
        if (I_Sta<0)
            continue     % Unknown station; skip this observation
        end
        
        % Compute 4-way range and partials using interpolated trajectories
        % of user and relay satellite
        [rho,drho_dx_User,drho_dx_TDRS] = Range_4W(Obs(i).Mjd_UTC,TrjData(1)...
                                          ,TrjData(I_TDRS+1),Sta(I_Sta+1));
        
        % Residual (observed - computed)
        res = Obs(i).Range_4w - rho;
        
        % Partials of modelled observations
        dzdX   = zeros(n_est,1);
        offset = 8*I_TDRS;
        for k=1:8
            dzdX(k)        = drho_dx_User(k);   % Partials w.r.t. user satellite
            dzdX(offset+k) = drho_dx_TDRS(k);   % Partials w.r.t. relay satellite
        end
        
        % Statistics
        n_obs      = n_obs + 1;
        Sum_res    = Sum_res + res;
        Sum_ressqr = Sum_ressqr + res*res;
        
        % Process data equation
        LSQAccumulate(n_est, dzdX, res, sigma_range);
        
        % Output
        [year,mon,day,hr,min,sec] = invjday(Obs(i).Mjd_UTC+2400000.5);
        fprintf(' %4d/%2.2d/%2.2d  %2.2d:%2.2d:%6.3f', year, mon, day, hr, min, sec);
        fprintf('%6d', Obs(i).StaNo);
        fprintf('%6d', Obs(i).TdrsNo);
        fprintf('%13.4f', Obs(i).Range_4w/1000);
        fprintf('%12.4f', rho/1000);
        fprintf('%10.2f', Obs(i).Range_4w-rho);
        fprintf('\n');
    end
    Aux(1).Mjd_UTC = Mjd0_UTC;
    Aux(2).Mjd_UTC = Mjd0_UTC;
    Aux(3).Mjd_UTC = Mjd0_UTC;
    
    % Solve least-squares system
    dX = LSQSolve(n_est, R, d);
    SigX = LSQStdDev(R, n_est);
    X = X + dX';
    
    % Update epoch states and parameters
    for i_sat=0:n_sat-1
        offset = 8*i_sat;
        y0(:,i_sat+1)   = X(offset+1:offset+6);
        Aux(i_sat+1).Cd = X(offset+7);
        Aux(i_sat+1).Cr = X(offset+8);
    end
    
    % Residuals statistics
    Mean = Sum_res/n_obs;
    RMS  = sqrt(Sum_ressqr/n_obs-Mean*Mean);
    
    % Print results
    fprintf('%70s',' ');
    fprintf('_______\n');
    fprintf('%70s',' Mean ');
    fprintf('%.2f', Mean);
    fprintf(' m\n');
    fprintf('%70s',' RMS  ');
    fprintf('%.2f', RMS);
    fprintf(' m\n\n\n');
    fprintf('Results of iteration ');
    fprintf('%d', iterat);
    fprintf('\n\n');
    fprintf('  Parameter            ');
    fprintf('a priori   tot.corr. last corr.       final       sigma\n');
    
    for ii=1:8
        fprintf('  ');
        if (ii<4)
            fprintf('%s', Label(ii));
            fprintf(' [m]  ');
        elseif (3<ii && ii<7)
            fprintf('v');
            fprintf('%s', Label(ii-3));
            fprintf('[m/s]');
        elseif(6<ii && ii<9)
            fprintf('C');
            fprintf('%s', Label(ii-3));
            fprintf('     ');
        end
        fprintf(' User  ');
        fprintf('%14.2f', X_apr(ii));
        fprintf('%11.2f', X(ii)-X_apr(ii));
        fprintf('%11.2f', dX(ii));
        fprintf('%14.2f', X(ii));
        fprintf('%11.2f', SigX(ii));
        fprintf('\n');
    end
    
    for ii=9:16
        fprintf('  ');
        if (ii<12)
            fprintf('%s', Label(ii-8));
            fprintf(' [m]  ');
        elseif (11<ii && ii<15)
            fprintf('v');
            fprintf('%s', Label(ii-11));
            fprintf('[m/s]');
        elseif(14<ii && ii<17)
            fprintf('C');
            fprintf('%s', Label(ii-11));
            fprintf('     ');
        end
        fprintf(' TDRS-');
        fprintf('%d', TdrsNo(ceil(ii/8)));
        fprintf('%14.2f', X_apr(ii));
        fprintf('%11.2f', X(ii)-X_apr(ii));
        fprintf('%11.2f', dX(ii));
        fprintf('%14.2f', X(ii));
        fprintf('%11.2f', SigX(ii));
        fprintf('\n');
    end
    
    for ii=17:24
        fprintf('  ');
        if (ii<20)
            fprintf('%s', Label(ii-16));
            fprintf(' [m]  ');
        elseif (19<ii && ii<23)
            fprintf('v');
            fprintf('%s', Label(ii-19));
            fprintf('[m/s]');
        elseif(22<ii && ii<25)
            fprintf('C');
            fprintf('%s', Label(ii-19));
            fprintf('     ');
        end
        fprintf(' TDRS-');
        fprintf('%d', TdrsNo(ceil(ii/8)));
        fprintf('%14.2f', X_apr(ii));
        fprintf('%11.2f', X(ii)-X_apr(ii));
        fprintf('%11.2f', dX(ii));
        fprintf('%14.2f', X(ii));
        fprintf('%11.2f', SigX(ii));
        fprintf('\n');
    end
    fprintf('\n\n');
end

profile viewer
profile off

