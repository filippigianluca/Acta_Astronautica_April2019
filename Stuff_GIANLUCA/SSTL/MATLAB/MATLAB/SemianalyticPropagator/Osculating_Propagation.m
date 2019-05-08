%% ========================================================================
% FABLE: Osculating Propagator.
% Perturbations:
% - J2, J3, J4, J5
% - Atmospheric drag
% - Constant acceleration in rth reference frame
% - Constant tangential acceleration
% - Constant Inertial acceleration
% - Solar radiation pressure
% - Third body
% =========================================================================
% Marilena Di Carlo, 2015



clear 
clc
close all

addpath('ChebyshevCoefficients')
addpath('IntegralsJ3')
addpath(genpath('../'))


%% Constants

% Earth Gravitational Constant [DU^3/TU^2]
constants.mu = 1;

% Distance unit DU - Earth Radius [km]
constants.DU = 6378.136;

% Time unit TU [s]
constants.TU = 806.78;

% Gravity acceleration [m/s2]
constants.g0 = 9.8;

% Adimensional Earth Radius
constants.R = 1;



%% Input

% Plots of the solutions?
i_plot = 1;

% Initial orbital elements [DU and rad]
kep0 = [(7000)/constants.DU ...
        0.001 ...
        80*pi/180 ...
        10*pi/180 ...
        60*pi/180 ...
        200*pi/180];

% Final true longitude
L_end = 5.580234237;


% Perturbations?

% Atmospheric drag
flag_drag = 0;

% J2
flag_J2   = 1;

% J3
flag_J3 = 1;

% J4
flag_J4 = 0;

% J5
flag_J5 = 0;

% Constant acceleration in rth reference frame
flag_rth  = 0;

% Constant tangential acceleraiton
flag_t    = 0;

% Constant inertial acceleration
flag_In   = 0;

% Sun and Moon
flag_Sun = 0;
flag_Moon = 0;

% SRP
% flag_SRP = 0;


%% Chebyshev interpolation of the atmospheric density
% 
% % Expansion order?
% order = 4;
% 
% if order == 4
%     % =========================================================================
%     % Expansion of order 4
%     % =========================================================================
%     
%     % Vector to define the height limits for each region
%     % h_crossing = [150 250 350 500 1000 2000 3000 4000 5000 6000 7000 8000 384000];
%     h_crossing = [150 250 350 500 700 1000 2000 3000 4000 384000];
%     h_crossing = h_crossing / DU;
%     
%     load Chebyshev_5_150_250.mat
% %     load Newton_5_150_250.mat
%     coeff(:,1) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_250_350.mat
% %      load Newton_5_250_350.mat
%     coeff(:,2) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_350_500.mat
% %     load Newton_5_350_500.mat
%     coeff(:,3) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_500_700.mat
% %     load Newton_5_500_700.mat
%     coeff(:,4) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_700_1000.mat
% %     load Newton_5_700_1000.mat
%     coeff(:,5) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_1000_2000.mat
% %     load Newton_5_1000_2000.mat
%     coeff(:,6) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_2000_3000.mat
%     %     load Newton_5_1000_2000.mat
%     coeff(:,7) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_5_3000_4000.mat
%     %     load Newton_5_1000_2000.mat
%     coeff(:,8) = coefficients1;
%     clear coefficients1
%     
%     % At higher altitude: no drag!
%     coeff(:,9) = zeros(5,1);
%     
% elseif order == 5
%     % =========================================================================
%     % Expansion of order 5
%     % =========================================================================
%     
%     % Vector to define the height limits for each region
%     % h_crossing = [150 250 350 500 1000 2000 3000 4000 5000 6000 7000 8000 384000];
%     h_crossing = [150 250 350 500 700 1000 2000 3000 4000];
%     h_crossing = h_crossing / DU;
%     
%     load Chebyshev_6_150_250.mat
%     coeff(:,1) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_6_250_350.mat
%     coeff(:,2) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_6_350_500.mat
%     coeff(:,3) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_6_500_700.mat
%     coeff(:,4) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_6_700_1000.mat
%     coeff(:,5) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_6_1000_2000.mat
%     coeff(:,6) = coefficients1;
%     clear coefficients1
%     
% elseif order == 7
%      % =========================================================================
%     % Expansion of order 7
%     % =========================================================================
%     
%     % Vector to define the height limits for each region
%     % h_crossing = [150 250 350 500 1000 2000 3000 4000 5000 6000 7000 8000 384000];
%     h_crossing = [150 250 350 500 700 1000 2000 3000 4000];
%     h_crossing = h_crossing / DU;
%     
%     load Chebyshev_8_150_250.mat
%     coeff(:,1) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_8_250_350.mat
%     coeff(:,2) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_8_350_500.mat
%     coeff(:,3) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_8_500_700.mat
%     coeff(:,4) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_8_700_1000.mat
%     coeff(:,5) = coefficients1;
%     clear coefficients1
%     
%     load Chebyshev_8_1000_2000.mat
%     coeff(:,6) = coefficients1;
%     clear coefficients1
%     
% end
% %% ________________________________________________________________________
% =========================================================================
% INPUT START
% _________________________________________________________________________
% =========================================================================



% Initial Keplerian Elements
Equin0 = kep2eq(kep0);

L0 = Equin0(6);



%% Propagation length

% Define number of revolutions
rev = 1;

% Define number of points of solution for each revolution
n_points = 1;

%% Perturbations included 

% =========================================================================
% Tangential thrust
% =========================================================================
if flag_t 
    
    % Acceleration [m/s2]
    a_t = 1e-4;
    
    % Adimensional acceleration
    epsilon_t = a_t * 10^(-3) * ( TU^2  / DU);
    
    % Elevation angle
    beta_t = 0;

else
    epsilon_t = 0;
    beta_t    = 0;
end

% =========================================================================
% Constant Acceleration in rth reference frame
% =========================================================================

if flag_rth
    
    % Acceleration [m/s2]
    a_rth = 1e-4;
    
    % Adimensional acceleration
    epsilon_rth = a_rth * 10^(-3) * ( TU^2  / DU);
    
    % Elevation angle [rad]
    beta_rth = pi/6;
    
    % Azimuth angle [rad]
    alpha_rth = pi/2;


else
    epsilon_rth = 0;
    beta_rth = 0;
    alpha_rth = 0;
end

% =========================================================================
% Inertial Acceleration
% =========================================================================

if flag_In
    
    % Acceleration [m/s^2]
    a_Inert = 1e-4;
    
    % Adimensional acceleration
    epsilon_In = a_Inert * 10^(-3) * ( TU^2  / DU);
        
    % Elevation angle
    beta_In = 0;
    
    % Initial azimuth angle
    alpha_In = 0;  
    
else
    
    epsilon_In = 0;
    beta_In = 0;
    alpha_In = 0;
    
end

% =========================================================================
% J2, J3, J4, J5
% =========================================================================

if flag_J2
    J2 = 1082.63e-6;
else
    J2 = 0;
end


if flag_J3
    J3 = -2.5327e-6;
else
    J3 = 0;
end

if flag_J4
    J4 =  -1.6196e-6;
else
    J4 = 0;
end

if flag_J5
    J5 = -2.2730e-7;
else
    J5 = 0;
end

geopotential.J2 = J2;
geopotential.J3 = J3;
geopotential.J4 = J4;
geopotential.J5 = J5;


% =========================================================================
% Sun, Moon
% =========================================================================
if flag_Sun
else
    third_body.flag_Sun = 0;
end

if flag_Moon
else
    third_body.flag_Moon = 0;
end

% =========================================================================
% Drag - check this
% =========================================================================

if flag_drag 
    
    % Compute drag integral analytically
    drag.num_an = 1;
    drag.CD = 2.2;
    drag.A_m =  1e-8 / DU^2; %[DU^2/kg]
    drag.rho = 1;
    drag.coeff = coeff;
    drag.h_crossing = h_crossing;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Atmospheric density over the considered orbit
    
    % Initialize true anomaly in the range 0 to 2*pi
    theta = linspace(0, 2*pi, 100);
    
    % Semiparameter [DU]
    p = kep0(1) * (1-kep0(2)^2);
    
    % Height [DU] for each value of the true anomaly
    h = p ./ (1 + kep0(2) * cos(theta)) - R;
    
    % Apogee and perigee height [DU]
    h_perigee = kep0(1) * (1-kep0(2)) - R;
    h_apogee  = kep0(1) * (1+kep0(2)) - R;
    
    if h_perigee < h_crossing(1)
        warning('Can not handle orbit with perigee lower than 150 km')
        warning('Change initial orbital elements')
        return
    end
    
    % Identify region of the Chebyshev fitting which are crossed by the orbit
    [~,index1] = find(h_crossing < h_perigee);
    [~,index2] = find(h_crossing > h_apogee);
    
    % Define vector h_crossing for the current orbit
    % 0 element vector means that the orbit lies entirely inside one region
    % 1 element vector means that the orbit lies in two regions
    % 2 elements vector means that the orbit lies in three regions
    % .... and so on
    if ~isempty(index1) && ~isempty(index2)
        h_crossing_orbit = h_crossing(:, index1(end) + 1 : index2 - 1);
    elseif isempty(index2)
        h_crossing_orbit = h_crossing(:, index1(end) + 1 : end);
    end
    
    % Number of crossing of regions
    number_crossing = length(h_crossing_orbit);
    
    % Theta of region crossing
    theta_crossing = acos((1/kep0(2)) .*( (p./(h_crossing_orbit + 1)) -1));
    
    % Define number of the regions which have to be considered
    regions_crossed = index1(end)  : index2 -1;
    
    
    
    rho_exact = zeros(1,length(theta));
    rho_Chebyshev = zeros(1,length(theta));
    
    
    index = 1;
    
    for i = 1 : length(theta)
        
        %     % Exact value of the atmospheric density [kg/DU^3]
        rho_exact(i) = exponential_atm_model_DUTU(h(i));
        
        % Moving from perigee to apogee
        if theta(i) < pi
            
            % Only if index is lower than number crossing:
            if index <= number_crossing
                if h(i) < h_crossing_orbit(index)
                else
                    index = index + 1;
                end
            end
            
            coeff_current =  coeff(:,regions_crossed(1,index));
            
            rho_Chebyshev(i) = 0;
            
            n = size(coeff,1);
            
            for k = 1 : n
                
                rho_Chebyshev(i) = rho_Chebyshev(i) + coeff_current(n-k+1) * (h(i)).^(k-1);
                
            end
            
        % Moving from apogee to perigee
        else
            
            if index >= 2
                if h(i) > h_crossing_orbit(index-1)
                else
                    index = index - 1;
                end
            end
            
            coeff_current =  coeff(:,regions_crossed(1,index));
            
            rho_Chebyshev(i) = 0;
            
%             n = length(coeff);
            
            for k = 1 : n
                
                rho_Chebyshev(i) = rho_Chebyshev(i) + coeff_current(n-k+1) * (h(i)).^(k-1);
                
            end
         
        end
              
    end

    if i_plot == 1
        figure
        subplot(2,1,1)
        semilogy(theta,rho_exact,'b','LineWidth',2)
        hold on
        plot(theta,rho_Chebyshev,'r','LineWidth',1)
        xlabel('True anomaly [rad]')
        ylabel('Density [kg/DU^3]')
        legend('Real density','Chebyshev density')
        subplot(2,1,2)
        plot(theta,h*DU,'k','LineWidth',2)
        hold on
        xlabel('True anomaly [rad]')
        ylabel('Height [km]')
        for i = 1 : length(h_crossing)
            line([0 2*pi],[h_crossing(i)*DU h_crossing(i)*DU])
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    drag.CD = 0;
end




%% Settings for numerical propagation

drag_NUM = drag;
% drag_NUM.A_m = 1e-8 / DU^2;

% Inputs to the numerical integration function:

% Distance unit DU
input.DU = constants.DU;

% Adimensional gravitational constant
input.mu = 1;

% Adimensional Earth radius
input.R = 1;

% J2
input.J2 = J2;


% J3
input.J3 = J3;

% Tangential acceleration
input.epsilon_t = epsilon_t;

% Constant acceleration in rth reference frame
input.epsilon_rth = epsilon_rth;
input.beta_rth    = beta_rth;
input.alpha_rth   = alpha_rth;

% Inertial acceleration
input.epsilon_In = epsilon_In;
input.alpha_In   = alpha_In;
input.beta_In    = beta_In;
input.gamma0_In  = alpha_In + Equin0(6);

input.a = kep0(1);
input.e = kep0(2);


%% Analytical Propagation

% Initialize results  of analytical propagator
EquinoctialElements = Equin0;
time = 0;


tic
% ---------------------------------------------------------------------
% Analytical integration
% ---------------------------------------------------------------------
Equin0_AN = Equin0;
Equin0_NUM = Equin0';
t0 = 0;

L = linspace(Equin0_AN(6), L_end, n_points);

% [Equin,t] = AnEquin_all_forward_tang_m(L, Equin0, epsilon_t, beta_t, epsilon_rth, alpha_rth, beta_rth, epsilon_In, alpha_In, beta_In, muadim, J2, R, ill_flag)
[Equin,t] = AnEquin_all_forward_tang_m(L, Equin0_AN, epsilon_t, beta_t,...
    epsilon_rth ,alpha_rth ,beta_rth ,...
    epsilon_In ,alpha_In ,beta_In ,...
    constants.mu, geopotential, third_body, 1, drag, 0, 0);
toc

time = [t0 t];
EquinoctialElements = [Equin0 Equin];

%% Numerical Propagation

tic
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
% [L,EquinNUM]=ode113(@(L,x)equations_motion(L,x,input,drag_NUM),linspace(Equin0_NUM(6), Equin0_NUM(6) + 2*pi*rev, n_points*rev - rev + 1),[Equin0_NUM(1:5) 0],options);
[L,EquinNUM]=ode113(@(L,x)equations_motion(L,x,input,drag_NUM),...
    linspace(Equin0_NUM(6), L_end, n_points),...
[Equin0_NUM(1:5) 0],options);

toc

EquinoctialElementsNUM = [EquinNUM(:,1:5) L];
timeNUM = EquinNUM(:,6);


%% Results

% Absolute error equinoctial elements
delta     = EquinoctialElements - EquinoctialElementsNUM(:,1:6)';

% Relative error equinoctial elements
delta_rel = (EquinoctialElements - EquinoctialElementsNUM(:,1:6)') ./ EquinoctialElementsNUM(:,1:6)';

% Relative error semimajor axis w.r.t. initial semimajor axis
delta_rel_a0 = (EquinoctialElements(1,:) - EquinoctialElementsNUM(:,1)') ./EquinoctialElements(1,1);

% Absolute error time [TU]
delta_time = time - timeNUM';

% Relative error time [TU]
delta_time_rel = (time - timeNUM') ./ timeNUM';


% -------------------------------------------------------------------------
% Keplerian and cartesian elements
% -------------------------------------------------------------------------

keplerian_elements     = zeros(6,size(EquinoctialElements,2));
keplerian_elements_NUM = zeros(6,size(EquinoctialElements,2));

cartesian_elements     = zeros(6,size(EquinoctialElements,2));
cartesian_elements_NUM = zeros(6,size(EquinoctialElements,2));


h_a = zeros(size(EquinoctialElements,2),1);
h_p = zeros(size(EquinoctialElements,2),1);

h_a_NUM = zeros(size(EquinoctialElements,2),1);
h_p_NUM = zeros(size(EquinoctialElements,2),1);

for k = 1 : size(EquinoctialElements,2)
    
    keplerian_elements(:,k) = eq2kep(EquinoctialElements(1:6,k));
    keplerian_elements_NUM(:,k) = eq2kep(EquinoctialElementsNUM(k,:))';
    
    cartesian_elements(:,k) = kep2cart(keplerian_elements(:,k),constants.mu);
    cartesian_elements_NUM(:,k) = kep2cart(keplerian_elements_NUM(:,k),constants.mu);
    
    h_a(k) = keplerian_elements(1,k) * (1 + keplerian_elements(2,k)) - constants.R;
    h_p(k) = keplerian_elements(1,k) * (1 - keplerian_elements(2,k)) - constants.R;
    
    h_a_NUM(k) = keplerian_elements_NUM(1,k) * (1 + keplerian_elements_NUM(2,k)) - constants.R;
    h_p_NUM(k) = keplerian_elements_NUM(1,k) * (1 - keplerian_elements_NUM(2,k)) - constants.R;
    
    
end

height = keplerian_elements(1,:) .* (1-keplerian_elements(2,:).^2) ./...
    (1+keplerian_elements(2,:) .* cos(keplerian_elements(6,:))) - constants.R;
height_NUM = keplerian_elements_NUM(1,:) .* (1-keplerian_elements_NUM(2,:).^2) ./ ...
    (1+keplerian_elements_NUM(2,:) .* cos(keplerian_elements_NUM(6,:)))  - constants.R;

% Keplerian elements difference
delta_kep =( keplerian_elements - keplerian_elements_NUM) ;
delta_kep_rel =( keplerian_elements - keplerian_elements_NUM) ./ keplerian_elements_NUM;

% Cartesian elements difference
delta_cart = cartesian_elements - cartesian_elements_NUM;

%% Plots

if i_plot == 1
    
%     % ---------------------------------------------------------------------
%     % 4 plot figure: absolute error of (a,P1,Q1,T)
%     % ---------------------------------------------------------------------
%     figure
%     subplot(2,2,1)
% %     plot(linspace(0,rev*2*pi,size(EquinoctialElements,2)), abs(delta(1,:))*DU)
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), abs(delta(1,:))*constants.DU)
%     xlabel('L [deg]')
%     ylabel('\Delta a [km]')
%     grid on
%     subplot(2,2,2)
% %         plot(linspace(0,rev*2*pi,size(EquinoctialElements,2)), abs(delta(2,:)))
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), abs(delta(2,:)))
%     xlabel('L [deg]')
%     ylabel('\Delta P_1')
%     grid on
%     subplot(2,2,3)
% %         plot(linspace(0,rev*2*pi,size(EquinoctialElements,2)), abs(delta(4,:)))
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), abs(delta(4,:)))
%     xlabel('L [deg]')
%     ylabel('\Delta Q_1')
%     grid on
%     subplot(2,2,4)
% %         plot(linspace(0,rev*2*pi,size(EquinoctialElements,2)), abs(delta_time * TU ))
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), abs(delta_time * constants.TU ))
%     xlabel('L [deg]')
%     ylabel('\Delta T [s]')
%     grid on
%     
%     % ---------------------------------------------------------------------
%     % Absolute error of a
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), (delta(1,:))*constants.DU)
%     xlabel('L [deg]')
%     ylabel('\Delta a [km]')
%     grid on
% 
%     % ---------------------------------------------------------------------
%     % Absolute error of P1
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), delta(2,:))
%     xlabel('L [deg]')
%     ylabel('\Delta P_1')
%     grid on
% 
%     % ---------------------------------------------------------------------
%     % Absolute error of P2
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), delta(3,:))
%     xlabel('L [deg]')
%     ylabel('\Delta P_2')
%     grid on
% 
%     % ---------------------------------------------------------------------
%     % Absolute error of Q1
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), delta(4,:))
%     xlabel('L [deg]')
%     ylabel('\Delta Q_1')
%     grid on
% 
%     % ---------------------------------------------------------------------
%     % Absolute error of Q2
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), delta(5,:))
%     xlabel('L [deg]')
%     ylabel('\Delta Q_2')
%     grid on
% 
%     % ---------------------------------------------------------------------
%     % Absolute error time
%     % ---------------------------------------------------------------------
%     figure
%     plot(linspace(L0*180/pi,L_end*180/pi,n_points), delta_time * constants.TU )
%     xlabel('L [deg]')
%     ylabel('\Delta T [s]')
%     grid on



    % -------------------------------------------------------------------------
    % Semimajor axis and eccentricity variation
    % -------------------------------------------------------------------------
    figure
    subplot(1,2,1)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements(1,:)*constants.DU,'b','LineWidth',2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements_NUM(1,:)*constants.DU,'r','LineWidth',2)
    legend('Analytical Propagation','Numerical Propagation')
    xlabel('L [deg]')
    ylabel('a [km]')
    subplot(1,2,2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements(2,:),'b','LineWidth',2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements_NUM(2,:),'r','LineWidth',2)
    legend('Analytical Propagation','Numerical Propagation')
    xlabel('L [deg]')
    ylabel('e')
    
    
    % -------------------------------------------------------------------------
    % Angles
    % -------------------------------------------------------------------------
    figure
    subplot(1,3,1)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements(3,:)*180/pi,'b','LineWidth',2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements_NUM(3,:)*180/pi,'r','LineWidth',2)
    legend('Analytical Propagation','Numerical Propagation')
    xlabel('L [deg]')
    ylabel('i [deg]')
    subplot(1,3,2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements(4,:)*180/pi,'b','LineWidth',2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements_NUM(4,:)*180/pi,'r','LineWidth',2)
    legend('Analytical Propagation','Numerical Propagation')
    xlabel('L [deg]')
    ylabel('RAAN [deg]')
    subplot(1,3,3)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements(5,:)*180/pi,'b','LineWidth',2)
    hold on
    plot(linspace(L0*180/pi,L_end*180/pi,n_points),keplerian_elements_NUM(5,:)*180/pi,'r','LineWidth',2)
    legend('Analytical Propagation','Numerical Propagation')
    xlabel('L [deg]')
    ylabel('\omega [deg]')

end



%% 

% true_anomaly = linspace(0, 2*pi, 100);
% a = kep0(1);
% e = kep0(2);
% p = kep0(1) * (1-kep0(2)^2);
% r = p ./ (1+e*cos(true_anomaly));
% 
% 
% delta_a= J2 * constants.R^2 / kep0(1) * ( (kep0(1)./r).^3 - 1 /(1-e^2)^(3/2) + 1.5 .* sin(kep0(3))^2 .* ...
%     (-(kep0(1)./r).^3 +  1 /(1-e^2)^(3/2) + (kep0(1)./r).^3 .* cos(2*kep0(5)+2*true_anomaly) ) );
% 
% 
% delta_e = J2 * constants.R^2 / 4 * (-2/(kep0(1)^2 * e * sqrt(1-e^2)) + 2 * a * (1-e^2) ./ (e * r.^3) + sin(kep0(3))^2 * (3 /(a^2 * e * sqrt(1-e^2)) - 3 *a * (1-e^2) ./ (e * r.^3)+...
%                -3 * (1-e^2) * cos(true_anomaly + 2*kep0(5)) / p^2 - 3 * cos(2*true_anomaly + 2 * kep0(5)) / (a^2 * e * (1-e^2)) + ...
%                3*a*(1-e^2)*cos(2*(true_anomaly+kep0(5))) ./ (e*r.^3) - (1-e^2)*cos(3*true_anomaly + 2 * kep0(5))/p^2   ));
% 
% figure
% plot(true_anomaly, (kep0(1) + delta_a)*constants.DU,'LineWidth',2)
% hold on
% plot(true_anomaly, kep0(1)*constants.DU*ones(1,length(delta_a)),'k','LineWidth',2)
% xlabel('\theta [rad]')
% ylabel('a [km]')
% grid on
% 
% figure
% plot(true_anomaly, (kep0(2) + delta_e),'LineWidth',2)
% hold on
% plot(true_anomaly, kep0(2)*ones(1,length(delta_a)),'k','LineWidth',2)
% xlabel('\theta [rad]')
% ylabel('e []')
% grid on
% 

%%  Print final numerical and analytical equinoctial elements (for Romain C code validation)

% close all

p_num = keplerian_elements_NUM(1,end) * ...
    (1 - keplerian_elements_NUM(2,end)^2);
str_num = strcat('Numerical: ', num2str(p_num*constants.DU,14),',', ...
    num2str(EquinoctialElementsNUM(end,3),14), ',', ...
    num2str(EquinoctialElementsNUM(end,2),14), ',', ...
    num2str(EquinoctialElementsNUM(end,5),14), ',', ...
    num2str(EquinoctialElementsNUM(end,4),14), ',', ...
    num2str(EquinoctialElementsNUM(end,6),14));
disp(str_num)

p_an = keplerian_elements(1,end) * ...
    (1 - keplerian_elements(2,end)^2);
str_num = strcat('Analytical: ', num2str(p_an*constants.DU,14),',', ...
    num2str(EquinoctialElements(3,end),14), ',', ...
    num2str(EquinoctialElements(2,end),14), ',', ...
    num2str(EquinoctialElements(5,end),14), ',', ...
    num2str(EquinoctialElements(4,end),14), ',', ...
    num2str(EquinoctialElements(6,end),14));
disp(str_num)

