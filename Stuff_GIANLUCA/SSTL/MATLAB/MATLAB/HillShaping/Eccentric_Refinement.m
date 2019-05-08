clear all
clc
close all

% Add path to optimiser folder
addpath(genpath('MPAIDEA'))
addpath(genpath('Shaping'))
addpath(genpath('Results'))
addpath(genpath('Results\Eccentric'))
%% Load data and parameters

% Remember to update LB/UB if any changes

% Load the data
LoadData;
Opt = importdata('7700-8300_900-2000.dat');
t_dep = Opt.data(:,1);
TOF = Opt.data(:,2);
nr = Opt.data(:,3);
Lambda = Opt.data(:,6:11);

% step = [22 24 26 28]-1 ;
l = length(t_dep) ;  

%% Bodies selection

% Earth as departure body
DepBody = body(9);
ArrBody = body(10);


%% Optimization of trajectory

UB = [2e-1 1e-0 2e-1 1e-2 2*pi+1e-8 44.1*pi/180]; 
LB = [-2e-1 -1e-0 -2e-1 -1e-2 0 20*pi/180];


% Maximum number of function evaluations
nFeValMax = 5000;

%% Boundary conditions
lll = 0 ;
X = [];
DV = [];
INDICI = [];
T_viol = [];
for ii = 1:l
    
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
    departure_elements = Kep2Hill(Kep_Dep,mu_S);
    u_0 = departure_elements(4);
    % Obtain Hill Parameters from Kepler Parameters
    arrival_elements = Kep2Hill(Kep_Arr,mu_S);
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
    
    
    %             (arrival_elements(5:6)-departure_elements(5:6))*180/pi
    
    du_i_h = 2e-4;
    du_others = 2e-3;
    % Function to optimise - CHANGE THIS
    fitnessfcn = @(lambda) hill_refinement_DV(lambda,departure_elements,arrival_elements,TOF(ii),mu_S,shape_flag,param_flag,du_i_h,du_others,0);
%     keyboard
%     Lambda = rand(1,6) .* (UB - LB) + LB;
%     keyboard
    optionsFM = optimoptions('fmincon','Algorithm','interior-point','Display','Iter','TolX',1e-15,'TolFun',1e-6,'MaxFunEvals',2000,'MaxIter',250);
    [lambdaSol,f_val,exitFlag] = fmincon(fitnessfcn,Lambda(ii,:),[],[],[],[],LB,UB,[],optionsFM);
    
    X = [X ; lambdaSol];
    f_val_time = hill_shaping_DV(lambdaSol,departure_elements,arrival_elements,TOF(ii),mu_S,shape_flag,param_flag,du_i_h,du_others,1) ;
    DV = [DV ; f_val];
    T_viol = [T_viol ; -(f_val_time - f_val) ./ 0.04 ];
    INDICI = [INDICI ; [t_dep(ii) TOF(ii) nr(ii)]];
%     keyboard
    
end
    
Dep = INDICI(:,1); Tof = INDICI(:,2); Nr = INDICI(:,3); a2 = X(:,1) ; a4 = X(:,2) ; a6 = X(:,3) ;
G2 = X(:,4) ; G3 = X(:,5) ; Range = X(:,6) ;
T = table(Dep,Tof,Nr,DV,T_viol,a2,a4,a6,G2,G3,Range);
writetable(T,'Results\Eccentric\7700-8300_900-2000_refined.dat');
