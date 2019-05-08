% clear all; close all; clc;
addpath(genpath('DECOMPOSITION_ALGORITHM'));
addpath(genpath('MACSMINMAX_CURRENT')); 

addpath(genpath('sstl')); 
addpath(genpath('MATLAB'));


% addpath(genpath('/home/mmarchi/ESTECO/progetti/esa/iti/matlab_strathclyde/MINMAX_plus_DECOMPOSITION'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_decomposition:
%
% 1) function  'evaluate_minmax'
%                                    OR       load (minmax, minmin)
% 2) function  'evaluate_minmin'
%
% DECOMPOSITION:
%
% 3) function  'evaluate_Belief_coupled_vectors'
%
% 4) function  'Sampling_Belief_Plausibility'
%
% 5) function  'reconstruction_Belief_Plausibility'
%
% 6) function  'evaluate_Belief_Plausibility_exact'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset random numbers generator
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
%s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% INITIALISATION

init = str2func('init_algo_so');
init_algo_decomposition = str2func('init_algorithm');

%% Problem INITIALISATION
% HERE YOU CAN DEFINE:

% 1) THE TIPE OF OUTPUT (Belief, Plausibility or both);
% 2) THE TIPE OF INPUT  (Run minmax and minmin  OR load them);
% 3) IF YOU WANT TO EVALUATE THE EXACT CURVE;
% 4) THE NUMBER OF SUB-SYSTEMS IN WHICH THE FUNCTION IS DECOMPOSED IN;
% 5) THE NUMBER OF SAMPLES: the bigger is the number, the more accurate is the decomposition;
% 6) THE DIMENTIONS OF COUPLED AND UNCOUPLED VECTORS;
% 7) LOWER AND UPPER BOUNDS OF DESIGN VECTOR (d) AND UNCERTAIN VECTOR (u), AND DIMENTION OF d;
% 8) NUMBER OF OBJECT FUNCTIONS;
% 9) FUNCTION(s);
% 10)MAX NUMBER OF EVALUATIONS FOR THE OPTIMIZER. 



%% 1) choose the output
% in.output = 0 --> only Belief (from min-max)
% in.output = 1 --> only Plausibility (from min-min)
% in.output = 2 --> Belief and Plausibility

in.output = 0;

%% 2) choose the input
% in.input = 0 --> do minmax and minmin
% in.input = 1 --> load d, u_min, u_max
% in.input = 2 --> load d, run max and min

in.input = 0;

%% 3) do exact Belief
% in.exact_curve(s) = 0 --> no
% in.exact_curve(s) = 1 --> yes

in.exact_curves = 0;

%% 4)number of sub-functions decomposition

num_functions = 6;  % number of sub-functions in which the problem is decomposed 

%% 5) number of samples
num_samples = [1 1 1 1 1 5 5 5 5 1 1 1 1 1 1];    % number of samples for each Belief and Plausibility curve of coupled vector





in.num_functions = num_functions;                                          
for i = 1:in.num_functions/2*(in.num_functions-1)
    in.num_samples{i}=num_samples(i);
end

%% EPISTEMIC VECTOR u (focal elements)

in.dim_u = [0 0 0 0 0 1 ...  % u1, u2, u3, u4 ,u5 ,u6
              0 0 0 0 6 ...  % u12,u13,u14,u15,u16
                0 0 0 6 ...  % u23,u24,u25,u26            
                  0 0 6 ...  % u34,u35,u36
                    0 6 ...  % u45,u46
                      6 ...  % u56 
            ];    
                                            
% 7) bounds of uncertain vector u;  

% dimension (1) satellite, (2) subsystems, (3) mass subsystems, (4) Q_in,
% (5) Q_ex

% delta 
D_axis = 20;
D_e =    0.0012;
D_I =    0.07;
D_RAAN = 30;
D_peri = 0.5;
D_th =   0.025;
D_eff =  0.08;

% avarage values

axis = [68500.3 73250.2 86065.5 49646.4 42049];
e =    [0.902 0.77866 0.513 0.15392 0.001];
I =    [22.81 9.12 1.09 0.36 0.05];
RAAN = [86.63 86.79 85.96 86.85 270];
peri = [180.1 180.6 180.81 180.97 1];
th =   [0 180.08 180.84 4.25 359.95];
eff =  0.8;

in.lb_u{1} = {[eff-D_eff eff];... % efficiency
              [axis(1)-D_axis axis(1)]; [e(1)-D_e e(1)]; [I(1)-D_I I(1)]; [RAAN(1)-D_RAAN RAAN(1)]; [peri(1)-D_peri peri(1)]; [th(1)-D_th th(1)];... % orbit parameters fire 1
              [axis(2)-D_axis axis(2)]; [e(2)-D_e e(2)]; [I(2)-D_I I(2)]; [RAAN(2)-D_RAAN RAAN(2)]; [peri(2)-D_peri peri(2)]; [th(2)-D_th th(2)];... % orbit parameters fire 2
              [axis(3)-D_axis axis(3)]; [e(3)-D_e e(3)]; [I(3)-D_I I(3)]; [RAAN(3)-D_RAAN RAAN(3)]; [peri(3)-D_peri peri(3)]; [th(3)-D_th th(3)];... % orbit parameters fire 3
              [axis(4)-D_axis axis(4)]; [e(4)-D_e e(4)]; [I(4)-D_I I(4)]; [RAAN(4)-D_RAAN RAAN(4)]; [peri(4)-D_peri peri(4)]; [th(4)-D_th th(4)];... % orbit parameters fire 4
              [axis(5)-D_axis axis(5)]; [e(5)-D_e e(5)]; [I(5)-D_I I(5)]; [RAAN(5)-D_RAAN RAAN(5)]; [peri(5)-D_peri peri(5)]; [th(5)-D_th th(5)];... % orbit parameters fire 5
                         
              }; 

          
in.ub_u{1} = {[eff eff+D_eff];... % efficiency
              [axis(1) axis(1)+D_axis]; [e(1) e(1)+D_e]; [I(1) I(1)+D_I]; [RAAN(1) RAAN(1)+D_RAAN]; [peri(1) peri(1)+D_peri]; [th(1) th(1)+D_th];... % orbit parameters fire 1
              [axis(2) axis(2)+D_axis]; [e(2) e(2)+D_e]; [I(2) I(2)+D_I]; [RAAN(2) RAAN(2)+D_RAAN]; [peri(2) peri(2)+D_peri]; [th(2) th(2)+D_th];... % orbit parameters fire 2
              [axis(3) axis(3)+D_axis]; [e(3) e(3)+D_e]; [I(3) I(3)+D_I]; [RAAN(3) RAAN(3)+D_RAAN]; [peri(3) peri(3)+D_peri]; [th(3) th(3)+D_th];... % orbit parameters fire 3
              [axis(4) axis(4)+D_axis]; [e(4) e(4)+D_e]; [I(4) I(4)+D_I]; [RAAN(4) RAAN(4)+D_RAAN]; [peri(4) peri(4)+D_peri]; [th(4) th(4)+D_th];... % orbit parameters fire 4
              [axis(5) axis(5)+D_axis]; [e(5) e(5)+D_e]; [I(5) I(5)+D_I]; [RAAN(5) RAAN(5)+D_RAAN]; [peri(5) peri(5)+D_peri]; [th(5) th(5)+D_th];... % orbit parameters fire 5
                           
              };    
          
in.bpa{1} =  {[.5 .5];...   % efficiency
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...     
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...
              [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5]; [.5 .5];...        
              
              };   
          



%% DESIGN VECTOR d

in.lb_d = [35; 0];     % [D_V, cell] 5000
in.ub_d = [45.1; 1];   %             10000


dim_d = length(in.lb_d);           % 7) dimentions of the vector d, lower and upper bounds;
in.dim_d = dim_d;


%% FIXED PARAMETERS

problem.fix.day = 59;%[];
problem.fix.MJD0 = 6.939792e+03 -1; %58484+7/24-1;
problem.fix.DoD_max = 0.2;
problem.fix.SoC_max = 0.9;

in.fix = problem.fix;
problem.par_objfun{1}.fix = in.fix;
%% FUNCTION 

% type of problem
problem.sign_inner = 1;               % -1 will run minmin(NOT CHANGE)

% objectives
problem.n_obj = 1;                         % 8) number of objective functions;
problem.objfun =     {@sstl_fun};          % 9) objective function(s); 
problem.constraints ={[]};%{@sstl_constraints};       
% problem.par_objfun = {struct};

% design variables
                                                             
problem.dim_d = in.dim_d;                                           
problem.lb_d = in.lb_d;                                             
problem.ub_d = in.ub_d;                                             

% uncertain variables

problem.dim_u_i = in.dim_u;                                         
dim_u = sum(problem.dim_u_i);

problem.dim_u = dim_u;

for n =1:problem.n_obj

  
problem.lb_u{n} = in.lb_u{n}; 
problem.ub_u{n} = in.ub_u{n}; 


end



% End of problem initialisation

%% Algorithm parameter INITIALISATION

%% Decomposition algorithm parameters
[ algo_decomposition ] = init_algo_decomposition(problem);
algo_decomposition.par.nFeValMax = 4;
algo_decomposition.par.n_agents = 6;

%% Min-Max/Min-Min algorithm parameters
[ algo_minmax, algo_outer, algo_inner ] = init(problem);
%% Local search flags of super, outer and inner algorithms (set up by user through mF GUI)
algo_minmax.par_minmax.local_search_flags.validation = false;
algo_minmax.par_minmax.local_search_flags.inner = false;
algo_minmax.par_minmax.local_search_flags.outer = false;
%% Function evaluations of outer and inner algorithms (set up by user through mF GUI)
algo_outer.par.nFeValMax = 600;
algo_inner.par.nFeValMax = 240;
%% Number of agents of outer and inner algorithms (set up by user through mF GUI)
algo_outer.par.n_agents = 5;
algo_inner.par.n_agents = 6;

%% SCENARIOS

if in.input == 0                                     %  run minmax and minmin
   
    minmax = [];
    minmin = [];
    
    if in.output == 0 || in.output == 2
        
        % MIN-MAX
        [minmax] = evaluate_minmax(problem, algo_minmax, algo_outer, algo_inner);
        n_points_tot = length(minmax.d(:,1));   
    end
    
    if in.output == 1 || in.output == 2
        
        % MIN-MIN
        [minmin] = evaluate_minmin(problem, algo_minmax, algo_outer, algo_inner);
        n_points_tot = length(minmin.d(:,1));         
    end
    

    

    n_obj_tot = problem.n_obj;
    
    minmax_input = minmax;
    minmin_input = minmin;
    problem_input = problem;
    
elseif in.input == 1                                     % input d and u

    minmax = [];
    minmin = [];
    
    load('input_minmax_minmin');
    
    n_points_tot = length(minmax.d(:,1));
    n_obj_tot = problem.n_obj;
    
    minmax_input = minmax;
    minmin_input = minmin;
    problem_input = problem;
    
elseif in.input == 2                                     % input d, run u

    minmax = [];
    minmin = [];
    max = [];
    min = [];
    
    load('input_d.mat')  % d = [...]
    
    if in.output == 0 || in.output == 2
        [max] = evaluate_max(problem, d, algo_inner);
    end

    if in.output == 1 || in.output == 2
        [min] = evaluate_min(problem, d, algo_inner);
    end
    
    n_points_tot = 1;
    n_obj_tot = 1;
    
    minmax_input = max;
    minmin_input = min;
    problem_input = problem;
    
end

for n_point = 1:n_points_tot                  % all the design points
    
    for n_obj = 1:n_obj_tot
        
        if isfield(minmax_input,'d')
            minmax.d = minmax_input.d(n_point,:);
            minmax.u = minmax_input.u{n_obj}(n_point,:);
            minmax.f = minmax_input.f(n_point, n_obj);
        end
        if isfield(minmin_input,'d')
            minmin.d = minmin_input.d(n_point,:);
            minmin.u = minmin_input.u{n_obj}(n_point,:);
            minmin.f = minmin_input.f(n_point, n_obj);             
        end
        

        

        
        global nfevalglobal;
        nfevalglobal = 0;
        
        %%
        global num_maximization_decomposition;
        num_maximization_decomposition = 0;
        
        [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(in, problem, minmax, minmin, n_obj, n_point, algo_decomposition);
        
        %% sampling:
        
        [Sample] = Sampling_Belief_Plausibility(Partial_curve, decomposition, in);
        
        %% reconstruction
        
        [decomposition_end, Plot_decomposition, num_sample_tot, LIST] = reconstruction_Belief_Plausibility(in, problem, minmax, minmin, Sample, n_obj, n_point, algo_decomposition);
        
        savefolder = strcat('RESULTS/2subsystems/');
        mkdir(savefolder);
        save(strcat(savefolder,'simple_TC'));
        
        %% EXACT BELIEF CURVE: run one optimization for each focal element
        
        if in.exact_curves == 1
            
            global num_maximization_Belief_exact;
            num_maximization_Belief_exact = 0;
            [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(in, problem, minmax, minmin, n_obj, algo_decomposition);
            
            figure
            hold on
            
            
            if in.output == 0
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                stairs(LIST.F_Bel, LIST.Bel,'.','linewidth',2)
                legend('Belief Exact','Belief Decomposition')
            elseif in.output == 1
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
                stairs(LIST.F_Pl, LIST.Pl,'.','linewidth',2)
                legend('Plausibility Exact','Plausibility Decomposition')
            elseif in.output == 2
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
                
                stairs(LIST.F_Bel, LIST.Bel,'.','linewidth',2)
                stairs(LIST.F_Pl, LIST.Pl,'.','linewidth',2)
                legend('Belief e Plausibility Exact','Belief e Plausibility Decomposition')
            end
            legend('Belief exact','Decomposition')
            hold off
            
            savefolder = strcat('RESULTS/2subsystems_exact/');
            mkdir(savefolder);
            save(strcat(savefolder,'simple_TC_exact'));
            
        end
        
    end
    
end