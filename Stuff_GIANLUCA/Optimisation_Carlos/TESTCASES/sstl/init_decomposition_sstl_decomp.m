function [in, problem] = init_decomposition_sstl_decomp()


%% choose the output
% in.output = 0 --> only Belief
% in.output = 1 --> only Plausibility
% in.output = 2 --> Belief and Plausibility

in.output = 0;

%% choose the input
% in.input = 0 --> do minmax and minmin
% in.input = 1 --> load d, u_min, u_max

in.input = 1;

%% NUMBER OF SUBFUNCTIONS decomposition

num_functions = 2;
in.num_functions = num_functions;                                          % number of subfunctions in the decomposition

%% number of samples for each Belief and Plausibility curve of coupled vector

for i = 1:in.num_functions/2*(in.num_functions-1);
    in.num_samples{i}=3;
end


%% output of sensibilityAnalysis. it gives positions of [u(1) u(2) u(12)]

% SensitivityAnalysis_i(1) = [1];                 % not finished
% SensitivityAnalysis_i(2) = [2];
% SensitivityAnalysis_i(3) = [3];
% SensitivityAnalysis_i(4) = [4];
% SensitivityAnalysis_i(5) = [5];
% SensitivityAnalysis_i(6) = [6];
% 
% in.SensitivityAnalysis = [SensitivityAnalysis_i(1:6)];                     %   HC


%% EPISTEMIC VECTOR u (focal elements)

% in.dim_u = [3 2 5 0 3 3];    
% in.lb_u{1} = {[0.005, 0.01]; [0.034, 0.0885]; [0.5, 0.6];   [0.6, 0.8]; [0.1, 0.2]; [0.025, 0.03];    [0.8, 0.85]; [3.5, 3.6]; [0, 10];    [0, 10];  [0.5, 1];    [2, 2.2];    [-10, 5];  [1, 3];   [0.1, 0.5]; [0.5, 1.5];};
% in.ub_u{1} = {[0.009, 0.02]; [0.0884, 0.15]; [0.59, 0.7]; [0.79, 0.9];  [0.3, 0.5]; [0.0275, 0.0375]; [0.84, 0.9]; [3.59, 4];  [9.99, 20]; [9.99, 15]; [0.99, 1.5]; [2.19, 2.5]; [4.99, 10]; [2.9, 5]; [0.49, 1];  [1.49, 2];};
% in.bpa{1} = {[0.5 0.5]; [0.5 0.5]; [0.5 0.5]; [0.2 0.8]; [0.4 0.6]; [0.8 0.2]; [0.4 0.6]; [0.5 0.5]; [0.3 0.7]; [0.4 0.6]; [0.5 0.5]; [0.5 0.5]; [0.4 0.6]; [0.3 0.7]; [0.3 0.7]; [0.4 0.6];};
% in.dim_u = [22 10 6];  %  ONLY ENERGY(coupled: effincecy single, eff triple, assembly effic single, assembly effi triple)
% [ in.lb_u, in.ub_u ] = bound(  );
% in.lb_u{1}{33} = [16.42,19]./100; % effincecy single
% in.ub_u{1}{33} = [18.99,21.16]./100;
% in.lb_u{1}{34} = [21.7,25]./100;
% in.ub_u{1}{34} = [24.99,29.81]./100;
% in.lb_u{1}{35} = [0.96,0.98]; % assemply single
% in.ub_u{1}{35} = [0.98,0.99];
% in.lb_u{1}{36} = [0.9,0.92]; % assembly triple
% in.ub_u{1}{36} = [0.92,0.94];
% in.lb_u{1}{37} = [0.94,0.96]; % efficiency arrey-battery-load
% in.ub_u{1}{37} = [0.96,0.98];
% in.lb_u{1}{38} = [0.95,0.97]; % efficiency arrey-load
% in.ub_u{1}{38} = [0.97,0.99];
% in.bpa{1} = repmat({[0.5 0.5]},[34,1]);
%% DESIGN VECTOR d

in.dim_d = 2;
in.lb_d = [0.0 ; 0.150];
in.ub_d = [1.0 ; 0.255];
%%%%%%%

% type of problem
problem.sign_inner = 1; % -1 will run minmin
% objectives
problem.n_obj = 1;
problem.objfun = {@fsstl2_dp_decomp};

%% DESIGN VECTOR d



% different exchange
problem.maxnfeval = 5e4;


load ('orbits.mat')
info_u = [...
1.737704918    , .05;
1.737704918    , .05;
0.3278688525   , .02;
0.3278688525   , .02;
0.7106598985   , .05;
0.7106598985   , .05;
0.7106598985   , .05;
6.2944162437   , .02;
1.6393442623   , .02;
1.6393442623   , .02;
1.6393442623   , .02;
0.835          , .05;
0.835          , .05;

1.3            , .02;
1.3            , .02;
59.8984771574  , .05;
4.0983606557   , .02;
8.1218274112   , .02;
8.1218274112   , .02;
25.3807106599  , .10;
45.6852791878  , .10;
3.0456852792   , .10;
3.0456852792   , .10;
6.0913705584   , .10;
6.0913705584   , .10];

for i  =1:size(info_u,1)
    in.lb_u{1}{i,1} = [info_u(i,1), info_u(i,1)*(1+0.5*info_u(i,2))];
    in.ub_u{1}{i,1} = [info_u(i,1)*(1+0.5*info_u(i,2)), info_u(i,1)*(1+info_u(i,2))];

    if info_u(i,2) == 0.02
        in.lb_u{1}{i,1} = [info_u(i,1)];
        in.ub_u{1}{i,1} = [info_u(i,1)*(1+info_u(i,2))];    
    end

end
in.ub_u{1}{i+1,1} = [1.0-0.19,1.0-0.1642]; % efficiency single
in.lb_u{1}{i+1,1} = [1.0-0.2116,1.0-0.19];

in.ub_u{1}{i+2,1} = [1.0-0.25,1.0-0.217]; % efficiency triple
in.lb_u{1}{i+2,1} = [1.0-0.2981,1.0-0.25];

in.ub_u{1}{i+3,1} = [1.0-0.95,1.0-0.92]; % assembly single
in.lb_u{1}{i+3,1} = [1.0-0.99,1.0-0.95];

in.ub_u{1}{i+4,1} = [1.0-0.95,1.0-0.93]; % assembly triple
in.lb_u{1}{i+4,1} = [1.0-0.97,1.0-0.95];

in.ub_u{1}{i+5,1} = [1.0-0.96,1.0-0.94]; % efficiency night
in.lb_u{1}{i+5,1} = [1.0-0.98,1.0-0.96];

in.ub_u{1}{i+6,1} = [1.0-0.97,1.0-0.95]; % efficiency day
in.lb_u{1}{i+6,1} = [1.0-0.99,1.0-0.97];

in.bpa{1} = cell(size(in.lb_u{1},1),1);
for i = 1:size(in.lb_u{1},1)
    n_int = size(in.lb_u{1}{i,1},2);
    in.bpa{1}{i,1} = repmat(1/n_int,[1, n_int]);
end

groups = cell(1,20);
groups {1} = 1:4;
groups {2} = 5:7;
groups {3} = 8;
groups {4} = 9:11;
groups {7} = 12:13;

groups {8} = 14:15;
groups {11}= 16;
groups {13}= 17;
groups {15}= 18;
groups {16}= 19;
groups {17}= 20;
groups {19}= 21;
groups {20}= [22,24];
groups {21}= [23,25];

minimum = [6 6 6 6 6 6 6 6 6 6 6 6 6 6];
realistic = [1 19 15 1 5 19 19 1 19 15 1 5 19 19];
goal = [6 6 20 6 20 20 6 6 6 20 6 20 20 6];

par.day = goal;
par.orbits = orbits;
par.groups = groups;

problem.par_objfun = {par};

in.dim_u = [13 12 6];

% design variables
%dim_d = 2;                                                                % dim_d = 2;
problem.dim_d = in.dim_d;                                           % problem_minmax.dim_d = dim_d;
problem.lb_d = in.lb_d;                                             % problem_minmax.lb_d = repmat(-5,[dim_d,1]);
problem.ub_d = in.ub_d;                                             % problem_minmax.ub_d = repmat(5,[dim_d,1]);

% uncertain variables

problem.dim_u_i = in.dim_u;                                         % [1 1 1 1 1 1] from INPUT_DECOMPOSITION
dim_u = sum(problem.dim_u_i);

problem.dim_u = dim_u;

for n =1:problem.n_obj

% problem.lb_u{n} = repmat({in.lb_u{n}},[dim_u,1]);    
% problem.ub_u{n} = repmat({in.ub_u{n}},[dim_u,1]);   
problem.lb_u{n} = in.lb_u{n}; 
problem.ub_u{n} = in.ub_u{n}; 



% problem.lb_u = {repmat({in.lb_u_1},[dim_u,1]), repmat({in.lb_u_2},[dim_u,1])};   for multiobjective 
% problem.ub_u = {repmat({in.ub_u_1},[dim_u,1]), repmat({in.ub_u_2},[dim_u,1])};

end

 

end

    