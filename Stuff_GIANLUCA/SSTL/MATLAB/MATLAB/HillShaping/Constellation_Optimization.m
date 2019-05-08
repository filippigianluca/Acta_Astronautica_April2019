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


%% Bodies selection

% Earth as departure body
DepBody = body(6);

ArrBody = body(7);

% Departure time [MJD2000]
t_dep = [0];
l_dep = length(t_dep);

% Time of flight [days]
TOF = [500];
l_TOF = length(TOF);

%% Calculate feasible number of revolutions
n_max = TOF ./ (2.* pi ./ DepBody.n); % TOF / Period
n_min = TOF ./ (2.* pi ./ ArrBody.n); % TOF / Period

% keyboard
nr = linspace(n_min,n_max,2);
nr = 10 ; 
l_nr = length(nr); 





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
options.text = 1;



%% Lower and upper boundaries of the search space - CHANGE THIS

% The lower and upper boundaries of the search space depend on the given problem
% Constellation optimization
UB = [1e-8 1e-8 1e-8 1E3 2*pi 44*pi/180]; 
LB = [-1e-8 -1e-8 -1e-8 -1E-3 0 10*pi/180];

% Maximum number of function evaluations
nFeValMax = 1000;

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
% Initialize print file
Header1 = 'Dep' ; Header2 = 'TOF' ; Header3 = 'N' ; Header4 = 'DV' ; Header5 = 'T_viol';
Header6 = 'a2' ; Header7 = 'a4' ; Header8 = 'a6' ; Header9 = 'G3' ; Header10 = 'G4'; Header11 = 'Range';
fmt = '%d  %d  %d  %.4e  %.4e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e \n' ;
name = 'Results/Constellation.txt' ; 
fileId = fopen(name,'w');
fprintf(fileId, [ Header1 ' , ' Header2 ', ' Header3 '  ,  ' Header4 '   ,   ' Header5 '    ,   ' Header6 '    ,   ' Header7 '   ,   ' Header8 '    ,    ' Header9 '    ,    ' Header10 '     ,   ' Header11 '\n']);
fclose(fileId);
% keyboard

X = [];
MIN = [];
INDICI = [];
for iii = 1:l_dep
    for jjj = 1:l_TOF
        for kkk = 1:l_nr
            
            % True anomaly of Departure Body at Departure time
            u_dep = wrapTo2Pi(DepBody.n*(t_dep(iii)-DepBody.t0));
            u_dep = 0;
            Kep_Dep = [DepBody.a ; DepBody.e ; DepBody.i ; DepBody.Omega ; DepBody.om ; u_dep];
                        
            % True anomaly of Arrival body at (Departure Time + TOF)
            u_arr = wrapTo2Pi(ArrBody.n*(t_dep(iii)-DepBody.t0 + TOF(jjj)));
            u_arr = pi/4;
            Kep_Arr = [ArrBody.a ; ArrBody.e ; ArrBody.i ; ArrBody.Omega ; ArrBody.om ; u_arr];
           
            % Rotate around starting node-line to improve the accuracy
%             angle_rotation = pi/2-Kep_Dep(3);
%             
%             omega0 = Kep_Dep(4);
%             n1 = cos(omega0) ; n2 = sin(omega0) ; n3 = 0 ;
%             W = [0 -n3 n2 ; n3 0 -n1 ; -n2 n1 0];
%             Rot = eye(3) + sin(angle_rotation)*W+(2 * sin(angle_rotation/2)^2)*(W^2) ;
%             
%             Kep_Dep2 = Kep_Dep ;
%             Kep_Dep2(3) = Kep_Dep(3) + angle_rotation ;
%             
%             Kep_Arr2 = RF_forward_kep_rotation(Kep_Arr,omega0,angle_rotation);
%             keyboard
%             if max(abs(Kep_Arr2-Kep_Arr2))>1e-6
%                 error('Error in rotation')
%             end
            
            % Obtain Hill Parameters from Kepler Parameters
            departure_elements = Kep2Hill(Kep_Dep,mu_T);
            u_0 = departure_elements(4);
            % Obtain Hill Parameters from Kepler Parameters
            arrival_elements = Kep2Hill(Kep_Arr,mu_T);
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
            a0_g = 1/departure_elements(1);
            a1_g = (1/arrival_elements(1) - a0_g) / arrival_elements(4);
            lambda0 = [a0_g a1_g -1e-05 2.265670e-06 5.035454e+00];
%             lambda0 = zeros(1,5);
            [c,A,b] = Shaping_functions(departure_elements,arrival_elements,shape_flag,1,lambda0);
            x  = A\(b-c);
            [p] = Hill_Parameters_Ordering (x,lambda0,shape_flag,1);
            a0= p(1); a1= p(2); a2= p(3); a3 = p(4); a4 = p(5); a5 = p(6); a6=p(7);
            G0 = p(8); G1 = p(9); G2 = p(10); G3 = p(11);
            
            u_0 = departure_elements(4); u_f = arrival_elements(4) ;
            u = linspace(u_0,u_f,10000);
            r = 1./( a0 + a1.*u + a2.*u.^2 + ( a3 + a4.*u ).*cos(u) + ( a5 + a6.*u ).*sin(u) );
            figure
            plot(u,r,'r')
            grid on
            keyboard
% %             (arrival_elements(5:6)-departure_elements(5:6))*180/pi
            
            du_i_h = 3e-2;
            du_others = 5e-2;
            % Function to optimise - CHANGE THIS
            fitnessfcn = @(lambda) Constellation_hill_DV(lambda,departure_elements,arrival_elements,TOF(jjj),mu_S,shape_flag,param_flag,du_i_h,du_others,1);
          
            
%             t_out = tic ;
%             fitnessfcn([0 0 0 0 0 40*pi/180])
%             t_out = toc(t_out)
%             keyboard
%             
            % MP-AIDEA optimisation
            [x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options);
            [min_fval,ind] = min(fval) ; 
            MIN_step = Constellation_hill_DV(x(ind,:),departure_elements,arrival_elements,TOF(jjj),mu_S,shape_flag,param_flag,du_i_h,du_others,1) ; 
            T_viol = (min_fval - MIN_step) ./ 0.04 ;
            
            fileId = fopen(name,'a');
            fprintf(fileId, fmt , [t_dep(iii) TOF(jjj) nr(kkk) MIN_step T_viol x(ind,:)]);
            fclose(fileId);
            
        end
    end
    
end

t_out = tic ;
fitnessfcn(x(1,:))
t_out = toc(t_out)
