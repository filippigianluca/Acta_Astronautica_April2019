%Gtoc2_test
%12/10/15
%-------------------------------------------------------------------------%
% Apply the Spherical Shaping method to the GTOC2 competition, trying to
% reproduce the results
% For now, the algorithm works if at least the first leg is shaped; if it's
% not, an error message will appear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

%--------------------------------------------------------------------------
% Clean screen,memory and figures
%--------------------------------------------------------------------------
clc;
clear all;
close all;

%--------------------------------------------------------------------------
% Add folders containign subfunctions
%--------------------------------------------------------------------------
% Add folders containing subfunctions
%addpath(genpath('SphericalShapingFunctions'));
addpath('SphericalShapingFunctions');
addpath('TrajectoryPlotting');
addpath('OptimalityFunctions');

%--------------------------------------------------------------------------
% Global variables
%--------------------------------------------------------------------------
global d2r;
global r2d;
global AU2km;
global km2AU;
global T_sid;
global mu_S;
global year;

d2r = pi/180;
r2d = 180/pi;
AU2km = 1.49597870691e+08; 
km2AU = 1/AU2km;
year = 365.25;
T_sid = 86400;                                   % [s]
mu_S = 1.32712440018e11*km2AU^3*T_sid^2;         % [AU^3/day^2]

%% Earth and asteroids parameters

%--------------------------------------------------------------------------
%Earth (values from GTOC2 files where aavailable and JPL ephemeris)
%--------------------------------------------------------------------------
Earth.name = 'Earth';
Earth.a = 0.999988049532578;
Earth.e = 1.671681163160e-02;
Earth.i = 0.8854353079654e-03*d2r;
Earth.Omega = 175.40647696473*d2r;
Earth.w = 287.61577546182*d2r;
Earth.M0 = 257.60683707535*d2r;
Earth.n = sqrt(mu_S/Earth.a^3);     
Earth.t0 = 54000;

%--------------------------------------------------------------------------
% Set the parameters of the spacecraft
%--------------------------------------------------------------------------
m_i = 1500;                     % kg
Ceff = 4000*9.80665*1e-3;       % Ceff = Isp*g0 [km/s]

%--------------------------------------------------------------------------
% Load file with asteroids data (ATTENTION: ANGLE IN DEGREES)
%--------------------------------------------------------------------------
[data,f] = readFileImportdata('ast_ephem.txt');

%--------------------------------------------------------------------------
% Load database with asteroids combinations from GTOC2 solutions and select
% the desired combination
%--------------------------------------------------------------------------
AstComb;

%--------------------------------------------------------------------------
% Select the combination of asteroids (1-12)
%--------------------------------------------------------------------------
choice = menu('Select the asteroids combination','1','2','3','4','5','6','7','8','9','10','11','12');
switch choice
    case 1
       comb = 1; 
    case 2
       comb = 2; 
    case 3
       comb = 3; 
    case 4
       comb = 4; 
    case 5
       comb = 5; 
    case 6
       comb = 6; 
    case 7
       comb = 7; 
    case 8
       comb = 8; 
    case 9
       comb = 9; 
    case 10
       comb = 10; 
    case 11
       comb = 11; 
    case 12
       comb = 12; 
end

% Select the chosen asteroids combination of the user
AST = ASTEROIDS(comb,:);
depDate_vec = depDate_vector(comb,:);
arrDate_vec = arrDate_vector(comb,:);

%--------------------------------------------------------------------------
% Search for the desired asteroids in the database and insert in the
% structure ArrBody
%--------------------------------------------------------------------------
for i=1:4
    for j=1:max(size(data))-1
        if data(j,1)==AST(i)
            ArrBody(i).name = num2str(AST(i));
            ArrBody(i).a = data(j,2);
            ArrBody(i).e = data(j,3);
            ArrBody(i).i = data(j,4)*d2r;
            ArrBody(i).Omega = data(j,5)*d2r;
            ArrBody(i).w = data(j,6)*d2r;
            ArrBody(i).M0 = data(j,7)*d2r;
            ArrBody(i).n = sqrt(mu_S/ArrBody(i).a^3);
            ArrBody(i).t0 = data(j,8);
        end
    end
end

%% Variables of Spherical Shaping (some of these variables are meaningful
% only if the Spherical Shaping Function with the Newton cycle is used,
% otherwise they are not used in the Spherical Shaping With FSolve

%--------------------------------------------------------------------------
% Set vector with number of revolutions
%--------------------------------------------------------------------------
Nr_vect = [0:4];

%--------------------------------------------------------------------------
% Set coefficient for Newton-Rapshon and plotIter for plotting a2-TOFcomp
% with iteration
%--------------------------------------------------------------------------
a2 = 0;
h = 0.0001;
plotIter = 0;

%--------------------------------------------------------------------------
% Set tolerance,step and max iterations
%--------------------------------------------------------------------------
toll = 1e-3;                  % at least 1e-3 (1 day of error for TOF=100d)
thetaStep = 2*pi/150;
maxIter = 500;

%% Call the code to compute the trajectory
disp('-------------------------------------------------------------------');
disp('SPHERICAL SHAPING METHOD APPLIED TO REPRODUCE A GTOC2 SOLUTION');
disp('-------------------------------------------------------------------');
fprintf('The 4 asteroids of the combination n° %g are: \n',comb);
fprintf('%g \t %g \t %g \t %g \n',AST(1),AST(2),AST(3),AST(4));
disp('-------------------------------------------------------------------');

%-------------------------------------------------------------------------%
% Initialize matrices (also multidimensional) containing the results
%-------------------------------------------------------------------------%
Vinf = [3.5 0 0 0]*km2AU*T_sid;
alpha = 0*d2r;
beta = 0*d2r;
DVLegs = [];
NrLegs = [];
thetaLegs = [];
RLegs = [];
PhiLegs = [];
oneLegNotConv = 0;

% Number of grid points for each leg
NPointsInt = 100;
%-------------------------------------------------------------------------%
% Start cycle on number of legs
%-------------------------------------------------------------------------%
tic

for Nleg=1:4
    
    % Set reference values
    convergenceForOneNr = 0;
    DV_min = 1e10;
    TOF = arrDate_vec(Nleg) - depDate_vec(Nleg);
    
    % Select DepartureBody
    if Nleg==1
        DepBody = Earth;
    else
        DepBody = ArrBody(Nleg-1);
    end
    
    % Cycle on Nr
    for k=1:length(Nr_vect)
        
        % Choose which algorithm to use
        %[thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShaping(DepBody,ArrBody(Nleg),depDate_vec(Nleg),TOF,Nr_vect(k),toll,thetaStep,maxIter,a2,h,plotIter,Vinf(Nleg),alpha,beta);
        [thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime,Ucart,Vcart] = SphericalShapingWithFSolve(DepBody,ArrBody(Nleg),depDate_vec(Nleg),TOF,Nr_vect(k),toll,thetaStep,a2,Vinf(Nleg),alpha,beta);
        
        % Save variables if loop has converged and if the DV is smaller than
        % the one of the previous leg (if previous leg converged)
        if isLoopConverged==1
            if DV < DV_min
                convergenceForOneNr = 1;
                DV_min = DV;
                Nr_min = Nr_vect(k);
                R_min = R;
                Phi_min = Phi;
                theta_min = thetaVec;
                cont_min = cont;
                u_mat_min = u_mat;
                TPrime_min = TPrime;
                Ucart_min = Ucart;
                Vcart_min = Vcart;
            end
        end
    end
    
    % Error if leg not converged (for any Nr)
    if convergenceForOneNr == 0
        fprintf('--------------------------------------------------------------- \n');
        disp(['ERROR - The loop for the leg N ' num2str(Nleg) ' has not converged; no trajectory found with the desired parameters']);
        fprintf('The number of iterations done is: %g / %g   \n',cont,maxIter);
        fprintf('--------------------------------------------------------------- \n');
        oneLegNotConv = 1;
        
        % Saving of trajectory if converged (for at least one Nr)
    else
        fprintf('--------------------------------------------------------------- \n');
        disp(['CONVERGENCE - The loop for the leg N '  num2str(Nleg)  ' has converged; a trajectory has been found']);
        fprintf('With a number of revolutions equal to:  %g   \n', Nr_min);
        fprintf('The number of iterations done is: %g / %g   \n',cont_min,maxIter);
        fprintf('--------------------------------------------------------------- \n');
        
        % Save DV,Nr
        DVLegs(Nleg) = DV_min;
        NrLegs(Nleg) = Nr_min;
        
        %---------------------------------------------------------------
        % Interpolate to have the same length for the vectors for all the legs
        %---------------------------------------------------------------
        thetaMin = linspace(theta_min(1),theta_min(end),NPointsInt);
        thetaLegs(Nleg,:) = thetaMin;
        RLegs(Nleg,:) = spline(theta_min,R_min,thetaMin);
        PhiLegs(Nleg,:) = spline(theta_min,Phi_min,thetaMin);
        TPrimeLegs(Nleg,:) = spline(theta_min,TPrime_min,thetaMin);
        
        % Cycle for the matrices of U and V for interpolation in 1 dim
        for i=1:3
            UcartLegs(Nleg,i,:) = spline(theta_min,Ucart_min(i,:),thetaMin);
            VcartLegs(Nleg,i,:) = spline(theta_min,Vcart_min(i,:),thetaMin);
        end
        
        % Obtain the cartesian state
        Xt = RLegs(Nleg,:).*cos(PhiLegs(Nleg,:)).*cos(thetaLegs(Nleg,:));
        Yt = RLegs(Nleg,:).*cos(PhiLegs(Nleg,:)).*sin(thetaLegs(Nleg,:));
        Zt = RLegs(Nleg,:).*sin(PhiLegs(Nleg,:));
        
        % Store it in a matrix
        RMatLegs(Nleg,1,:) = Xt;
        RMatLegs(Nleg,2,:) = Yt;
        RMatLegs(Nleg,3,:) = Zt;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Previous version, where reference length is the one of the first leg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             lengthTheta(Nleg) = length(theta_min);
        %             thetaLegs(Nleg,:) = theta_min;
        %             RLegs(Nleg,:) = R_min;
        %             PhiLegs(Nleg,:) = Phi_min;
        %             TPrimeLegs(Nleg,:) = TPrime_min;
        %             UcartLegs(Nleg,:,:) = Ucart_min;
        %             VcartLegs(Nleg,:,:) = Vcart_min;
        %
        %             % Obtain the cartesian state
        %             Xt = R_min.*cos(Phi_min).*cos(theta_min);
        %             Yt = R_min.*cos(Phi_min).*sin(theta_min);
        %             Zt = R_min.*sin(Phi_min);
        %
        %             % Store it in a matrix
        %             RMatLegs(Nleg,1,:) = Xt;
        %             RMatLegs(Nleg,2,:) = Yt;
        %             RMatLegs(Nleg,3,:) = Zt;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Plot the control as a function of angle theta (not time!!)
        figure;
        subplot(3,1,1);
        plot(theta_min*r2d,u_mat_min(1,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
        ylabel('Ut [km/s^2]');
        subplot(3,1,2);
        plot(theta_min*r2d,u_mat_min(2,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
        ylabel('Un [km/s^2]');
        subplot(3,1,3);
        plot(theta_min*r2d,u_mat_min(3,:)*AU2km/T_sid^2,'b','LineWidth',2.0);
        ylabel('Uh [km/s^2]');
        xlabel('Theta [deg]');
        tit = strcat('Control acceleration for leg N', num2str(Nleg));
        suptitle(tit);
        set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold');
    end
    
end
t4L = toc;

%% Results

%--------------------------------------------------------------------------
% Check if all the legs have converged, and if so, calculate the total DV,
% the final mass, and plot the trajectory
%--------------------------------------------------------------------------
if oneLegNotConv == 1
    fprintf('--------------------------------------------------------------- \n');
    fprintf('ERROR - At least one of the four legs has not converged for COMBINATION N %g\n',comb);
    fprintf('--------------------------------------------------------------- \n');
    
else
    fprintf('--------------------------------------------------------------- \n');
    fprintf('OK - All the four legs have converged for COMBINATION N %g\n',comb);
    disp('-------------------------------------------------------------------');
    fprintf('The total computational time is: %f  min \n',t4L/60);
    fprintf('--------------------------------------------------------------- \n');
    
    % DV tot and mass calculation
    DVtot = sum(DVLegs);
    m_f = m_i*exp(-DVtot*AU2km/T_sid/Ceff);
    
    % Print data
    disp('-------------------------------------------------------------------');
    fprintf('The DeltaV required and the final mass are \n');
    fprintf('Performance index, J: %f kg/yr \n', m_f/((arrDate_vec(end) - depDate_vec(1))/year));
    fprintf('Final-Initial mass: %f  / %f kg \n',m_f,m_i);
    fprintf('DeltaV of 4 legs: %f - %f - %f - %f -   km/s \n',DVLegs(1)*AU2km/T_sid,DVLegs(2)*AU2km/T_sid,DVLegs(3)*AU2km/T_sid,DVLegs(4)*AU2km/T_sid);
    fprintf('DeltaV: %f  km/s \n',DVtot*AU2km/T_sid);
    fprintf('--------------------------------------------------------------- \n');
    
end

%--------------------------------------------------------------------------
% Plot trajectory
%--------------------------------------------------------------------------
plotTrajectoryMultBodiesSpherNoCont(Earth,ArrBody,depDate_vec,arrDate_vec,thetaLegs,RLegs,PhiLegs);
plotTrajectoryMultBodiesCartContMultTraj(Earth,ArrBody,depDate_vec,arrDate_vec,RMatLegs,UcartLegs);

%% Optimal Control Part (Verification and Translation into Optimal Trajectory)

% Execute only if all the trajectory has been shaped
if (oneLegNotConv == 0)
    % Initialize vector containing time
    timeLegs = [];
    
    % Cycle on number of legs
    for Nleg=1:4
        
        % Compute TOF
        TOF = arrDate_vec(Nleg) - depDate_vec(Nleg);
        
        % First obtain the time history vector
        timeLegs(Nleg,:) = timeEvolution(thetaLegs(Nleg,:),TPrimeLegs(Nleg,:),depDate_vec(Nleg),TOF,toll);
        
    end
    
    % Call verify optimality function
    verifyOptimalityGTOC2PhasesSeparated(Earth,ArrBody,RMatLegs,VcartLegs,UcartLegs,timeLegs,mu_S,DVLegs,depDate_vec,arrDate_vec);
    
end