clear all
clc
close all

% Add path to optimiser folder
addpath(genpath('MPAIDEA'))
addpath(genpath('Shaping'))
addpath(genpath('Results'))
%% Load data and parameters

% Load the data
LoadData;


% Departure time [MJD2000]
t_dep = [7900:20:8200];
% t_dep = 8000 ;
l_dep = length(t_dep);

% Time of flight [days]
TOF = [900:20:1100];
% TOF = 1600 ;
l_TOF = length(TOF);

Destination_Output = 'Results\Eccentric\7900-8200_900-1100_1.dat';
nr = [1] ;
% nr = 1 ;
l_nr = length(nr); 

%% Bodies selection

% Earth as departure body
DepBody = body(9);
ArrBody = body(10);


%% Optimization of trajectory


%% Dimension of the problem - CHANGE THIS

D = 6;


%% Set the parameters for MP-AIDEA

% -------------------------------------------------------------------------
% Number of populations
% -------------------------------------------------------------------------
pop_number = 4;

% -------------------------------------------------------------------------
% Number of individuals in each population
% -------------------------------------------------------------------------
NP = D;

% -------------------------------------------------------------------------
% Dimension of the bubble for the local restart. If empty, MP-AIDEA will
% adapt it
% -------------------------------------------------------------------------
options.delta_local = [];
% options.delta_local = 0.1;

% -------------------------------------------------------------------------
% Distance from cluster of global minima for global restart of the
% population
% -------------------------------------------------------------------------
options.delta_global = 0.1;

% -------------------------------------------------------------------------
% Threshold for contraction of the population
% -------------------------------------------------------------------------
options.rho = 0.2;

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% -------------------------------------------------------------------------
options.max_LR = [];
% options.max_LR = 5;

% -------------------------------------------------------------------------
% Choose the Differential Evolution (DE) strategies. 
% -------------------------------------------------------------------------
% The DE of MP-AIDEA uses two DE strategies, with probability
% defined by options.prob_DE_strategy (see later).
% DE/Rand, DE/CurrentToBest and DE/Best are well know DE strategies.
% Uncomment the following line for DE strategies DE/Rand and DE/CurrentToBest:
options.DE_strategy = 1;
% Uncomment the following line for DE strategies DE/Rand and DE/Best:
% options.DE_strategy = 2;


% -------------------------------------------------------------------------
% Probability of having DE strategy 1 rather than DE strategy 2 (DE
% strategies 1 and 2 have been defined in options.DE_strategy)
% -------------------------------------------------------------------------
% Example: if options.DE_strategy was set to 1, options.prob_DE_strategy
% defines the probability of using DE/Rand rather than DE/CurrentToBest
options.prob_DE_strategy = 0.5;


% -------------------------------------------------------------------------
% Parameter for the adaptation of CR and F
% -------------------------------------------------------------------------
% Value of CR (if empty, MP-AIDEA-ALR adapt will adapt it during the process)
% options.CR = 0.5;
options.CR = [];

% Value of F (if empty, MP-AIDEA will adapt it during the process)
% options.F = 0.5;
options.F = [];

% If options.CR and options.F are empty, define CRF for adaptation of CR
% and F:
options.dd_CRF = 3;


% -------------------------------------------------------------------------
% Display text during run?
% -------------------------------------------------------------------------
options.text = 0;



%% Lower and upper boundaries of the search space - CHANGE THIS

% The lower and upper boundaries of the search space depend on the given problem
% Mars optimization
% UB = [1e-1 1e-1 1e-1 1e-1 2*pi 44*pi/180]; 
% LB = [-1e-1 -1e-1 -1e-1 -1e-1 0 10*pi/180];
UB = [1e-0 1e-0 1e-0 1e-1 2*pi 44*pi/180]; 
LB = [-1e-0 -1e-0 -1e-0 -1e-1 0 10*pi/180];

% Maximum number of function evaluations
nFeValMax = 5000;

%% MP-AIDEA inputs

% Maximum number of function evaluations
options.nFeValMax = nFeValMax;

% Solutions are saved not only when nFeValMax has been reached but also at
% fraction of nFeValMax. Define this fraction in options.record:
options.record = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

% Initialise populations
population = zeros(NP,D,pop_number);

for s = 1 : pop_number
    
    pop = lhsdesign(NP,D,'criterion','maximin').*repmat(UB-LB,NP,1)+repmat(LB,NP,1);
    population(:,:,s) = pop;
end


options.population = population;


%% Optimisation

%% Boundary conditions


X = [];
DV = [];
INDICI = [];
T_viol = [];

for iii = 1:l_dep
    for jjj = 1:l_TOF
        for kkk = 1:l_nr
            
            % True anomaly of Departure Body at Departure time
            M_dep = DepBody.M0 + DepBody.n*(t_dep(iii)-DepBody.t0); % Extract mean anomali at t0
            Theta_dep = M2teta(M_dep,DepBody.e); % Convert to true anomaly at t0
            Kep_Dep = [DepBody.a ; DepBody.e ; DepBody.i ; DepBody.Omega ; DepBody.om ; Theta_dep];
            Cart_Dep = Kep2Cart(Kep_Dep,mu_S);
            
            
            % True anomaly of Arrival body at (Departure Time + TOF)
            M_arr = ArrBody.M0 + ArrBody.n*(t_dep(iii)+TOF(jjj)-ArrBody.t0);
            Theta_arr = M2teta(M_arr,ArrBody.e);
            Kep_Arr = [ArrBody.a ; ArrBody.e ; ArrBody.i ; ArrBody.Omega ; ArrBody.om ; Theta_arr];
            Cart_Arr = Kep2Cart(Kep_Arr,mu_S);
            
            % Rotate around starting node-line to improve the accuracy
            angle_rotation = pi/2-Kep_Dep(3); 
            
            omega0 = Kep_Dep(4);
            n1 = cos(omega0) ; n2 = sin(omega0) ; n3 = 0 ;
            W = [0 -n3 n2 ; n3 0 -n1 ; -n2 n1 0];
            Rot = eye(3) + sin(angle_rotation)*W+(2 * sin(angle_rotation/2)^2)*(W^2) ;        
            
            Cart_Dep2(1:3) = Rot * Cart_Dep(1:3);
            Cart_Dep2(4:6) = Rot * Cart_Dep(4:6);
            Cart_Arr2(1:3) = Rot * Cart_Arr(1:3);
            Cart_Arr2(4:6) = Rot * Cart_Arr(4:6);
            
            
            Kep_Dep2 = Kep_Dep ; 
            Kep_Dep2(3) = Kep_Dep(3) + angle_rotation ;
            
            Kep_Arr2 = RF_forward_kep_rotation(Kep_Arr,omega0,angle_rotation);
                      
            
            % Obtain Hill Parameters from Kepler Parameters
            departure_elements = Kep2Hill(Kep_Dep2,mu_S);
            u_0 = departure_elements(4);
            % Obtain Hill Parameters from Kepler Parameters
            arrival_elements = Kep2Hill(Kep_Arr2,mu_S);
            % Arrival argument of latitude
            u_fr = arrival_elements(4);
            
            while u_fr<u_0
                u_fr = u_fr + 2*pi;
            end
            
            du = u_fr - u_0 ;
            u_f = u_fr + 2*pi*nr(kkk);
            arrival_elements(4) = u_f;
            
            shape_flag = 1;
            param_flag = 0;
            
            while abs(arrival_elements(5)-departure_elements(5)) > pi
                
                if arrival_elements(5) > departure_elements(5)
                    arrival_elements(5) = arrival_elements(5) - 2*pi ;
                else
                    arrival_elements(5) = arrival_elements(5) + 2*pi ;
                end
            end
            while abs(arrival_elements(6)-departure_elements(6)) > pi
                
                if arrival_elements(6) > departure_elements(6)
                    arrival_elements(6) = arrival_elements(6) - 2*pi ;
                else
                    arrival_elements(6) = arrival_elements(6) + 2*pi ;
                end
                
            end
            
            du_i_h = 2e-4;
            du_others = 2e-3;
            
            % Function to optimise - CHANGE THIS
            fitnessfcn = @(lambda) hill_shaping_DV(lambda,departure_elements,arrival_elements,TOF(jjj),mu_S,shape_flag,param_flag,du_i_h,du_others,0);
        
%             a = importdata('7700-8300_900-2000.dat');
%             lambda0 = a.data(44,6:end);
%             t1 = tic ;
%             fitnessfcn(lambda0)
%             tf = toc(t1)
%             keyboard
            
            % MP-AIDEA optimisation
            [x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options);
            [min_fval,ind] = min(fval) ;
            
            %% Save Output for Table compilation
            X = [X ; x(ind,:)];
            MIN_step = hill_shaping_DV(x(ind,:),departure_elements,arrival_elements,TOF(jjj),mu_S,shape_flag,param_flag,du_i_h,du_others,1) ; 
            DV = [DV ; MIN_step];
            t_viol = (min_fval - MIN_step) ./ 0.04 ;
            T_viol = [T_viol ; t_viol ];
            INDICI = [INDICI ; [t_dep(iii) TOF(jjj) nr(kkk)]];
            
        end
        
    end
    
end


Dep = INDICI(:,1); Tof = INDICI(:,2); Nr = INDICI(:,3); a2 = X(:,1) ; a4 = X(:,2) ; a6 = X(:,3) ;
G2 = X(:,4) ; G3 = X(:,5) ; Range = X(:,6) ;
T = table(Dep,Tof,Nr,DV,T_viol,a2,a4,a6,G2,G3,Range);
writetable(T,Destination_Output);


clear all

% Add path to optimiser folder
addpath(genpath('MPAIDEA'))
addpath(genpath('Shaping'))
addpath(genpath('Results'))
%% Load data and parameters

% Remember to update LB/UB if any changes

% Load the data
LoadData;

Opt = importdata('Results\Eccentric\7900-8200_900-1100_1.dat');

t_dep = Opt.data(:,1);
TOF = Opt.data(:,2);
nr = Opt.data(:,3);
Lambda = Opt.data(:,6:11);

l = length(t_dep) ;  

%% Bodies selection

% Earth as departure body
DepBody = body(9);
ArrBody = body(10);


%% Optimization of trajectory

UB = [2e-1 1e-0 2e-1 1e-2 2*pi 44*pi/180]; 
LB = [-2e-1 -1e-0 -2e-1 -1e-2 0 20*pi/180];

% Maximum number of function evaluations
nFeValMax = 5000;

%% Boundary conditions
lll = 0 ;
X = [];
DV = [];
INDICI = [];
T_viol = [];
for ii = 1: l
    
    % True anomaly of Departure Body at Departure time
    M_dep = DepBody.M0 + DepBody.n*(t_dep(ii)-DepBody.t0); % Extract mean anomali at t0
    Theta_dep = M2teta(M_dep,DepBody.e); % Convert to true anomaly at t0
    Kep_Dep = [DepBody.a ; DepBody.e ; DepBody.i ; DepBody.Omega ; DepBody.om ; Theta_dep];
    Cart_Dep = Kep2Cart(Kep_Dep,mu_S);
    
    
    % True anomaly of Arrival body at (Departure Time + TOF)
    M_arr = ArrBody.M0 + ArrBody.n*(t_dep(ii)+TOF(ii)-ArrBody.t0);
    Theta_arr = M2teta(M_arr,ArrBody.e);
    Kep_Arr = [ArrBody.a ; ArrBody.e ; ArrBody.i ; ArrBody.Omega ; ArrBody.om ; Theta_arr];
    Cart_Arr = Kep2Cart(Kep_Arr,mu_S);
    
    % Rotate around starting node-line to improve the accuracy
    angle_rotation = pi/2-Kep_Dep(3);
    
    omega0 = Kep_Dep(4);
    n1 = cos(omega0) ; n2 = sin(omega0) ; n3 = 0 ;
    W = [0 -n3 n2 ; n3 0 -n1 ; -n2 n1 0];
    Rot = eye(3) + sin(angle_rotation)*W+(2 * sin(angle_rotation/2)^2)*(W^2) ;
    
    Cart_Dep2(1:3) = Rot * Cart_Dep(1:3);
    Cart_Dep2(4:6) = Rot * Cart_Dep(4:6);
    Cart_Arr2(1:3) = Rot * Cart_Arr(1:3);
    Cart_Arr2(4:6) = Rot * Cart_Arr(4:6);
    
    
    Kep_Dep2 = Kep_Dep ;
    Kep_Dep2(3) = Kep_Dep(3) + angle_rotation ;
    
    Kep_Arr3 = Cart2Kep(Cart_Arr2,mu_S) ;
    Kep_Arr2 = RF_forward_kep_rotation(Kep_Arr,omega0,angle_rotation);
    
    % Obtain Hill Parameters from Kepler Parameters
    departure_elements = Kep2Hill(Kep_Dep2,mu_S);
    u_0 = departure_elements(4);
    % Obtain Hill Parameters from Kepler Parameters
    arrival_elements = Kep2Hill(Kep_Arr2,mu_S);
    % Arrival argument of latitude
    u_fr = arrival_elements(4);
    
    while u_fr<u_0
        u_fr = u_fr + 2*pi;
    end
    
    du = u_fr - u_0 ;
    u_f = u_fr + 2*pi*nr(ii);
    arrival_elements(4) = u_f;
    
    shape_flag = 1;
    param_flag = 0;
    
    while abs(arrival_elements(5)-departure_elements(5)) > pi
        
        if arrival_elements(5) > departure_elements(5)
            arrival_elements(5) = arrival_elements(5) - 2*pi ;
        else
            arrival_elements(5) = arrival_elements(5) + 2*pi ;
        end
    end
    while abs(arrival_elements(6)-departure_elements(6)) > pi
        
        if arrival_elements(6) > departure_elements(6)
            arrival_elements(6) = arrival_elements(6) - 2*pi ;
        else
            arrival_elements(6) = arrival_elements(6) + 2*pi ;
        end
        
    end
    
    du_i_h = 2e-4;
    du_others = 2e-3;
    
    % Function to optimise - CHANGE THIS
    fitnessfcn = @(lambda) hill_refinement_DV(lambda,departure_elements,arrival_elements,TOF(ii),mu_S,shape_flag,param_flag,du_i_h,du_others,0);

    optionsFM = optimoptions('fmincon','Algorithm','interior-point','Display','Iter','TolX',1e-15,'TolFun',1e-6,'MaxFunEvals',5000,'MaxIter',450);
    [lambdaSol,f_val,exitFlag] = fmincon(fitnessfcn,Lambda(ii,:),[],[],[],[],LB,UB,[],optionsFM);
    
    X = [X ; lambdaSol];
    f_val_time = hill_shaping_DV(lambdaSol,departure_elements,arrival_elements,TOF(ii),mu_S,shape_flag,param_flag,du_i_h,du_others,1) ;
    DV = [DV ; f_val];
    T_viol = [T_viol ; -(f_val_time - f_val) ./ 0.04 ];
    INDICI = [INDICI ; [t_dep(ii) TOF(ii) nr(ii)]];
%     
    
end
    
Dep = INDICI(:,1); Tof = INDICI(:,2); Nr = INDICI(:,3); a2 = X(:,1) ; a4 = X(:,2) ; a6 = X(:,3) ;
G2 = X(:,4) ; G3 = X(:,5) ; Range = X(:,6) ;
T = table(Dep,Tof,Nr,DV,T_viol,a2,a4,a6,G2,G3,Range);
writetable(T,'Results\Eccentric\7900-8200_900-1100_1_revised.dat');