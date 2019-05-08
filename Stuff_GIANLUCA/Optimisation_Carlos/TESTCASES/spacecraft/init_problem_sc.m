function [problem] = init_problem_sc()

% type of problem
problem.sign_inner = 1; % -1 will run minmin
% objectives
problem.n_obj = 1;
problem.objfun = {@Spacecraft};
problem.par_objfun = {struct};

%% DESIGN d

problem.dim_d = 10;
problem.lb_d = [10;30;7;0;0;0.15;135;3;1;0];
problem.ub_d = [60;90;10;1;1;0.3;145;5;3;1];

% UNCERTAINTY u
problem.dim_u = 16; 
problem.lb_u{1} = {[0.005, 0.01]; [0.034, 0.0885]; [0.5, 0.6]; [0.5, 1]; [2, 2.2]; [-10, 5]; ... % AOCS (l, A, q, m, CD, dI)
                   [0.6, 0.8]; [1, 3]; [0.1, 0.5]; [0.5, 1.5]; [0.1, 0.2]; ...                   % TTC  (eta, Gt, Lt, Lother, Mrfdn)
                   [0.025, 0.03]; [0.8, 0.85]; [3.5, 3.6]; [0, 10]; [0, 10] };                   % EPS  (Dcell, eta, rho, dP, Tmax)

problem.ub_u{1} = {[0.01, 0.02]; [0.0885, 0.15]; [0.6, 0.7]; [1, 1.5]; [2.2, 2.5]; [5, 10]; ...   % AOCS
                   [0.8, 0.9]; [3, 5]; [0.5, 1];  [1.5, 2]; [0.3, 0.5]; ...                      % TTC
                   [0.0275, 0.0375]; [0.85, 0.9]; [3.6, 4];  [10, 20]; [10, 15] };               % EPS
  
problem.bpa{1}  = {[0.5 0.5]; [0.5 0.5]; [0.5 0.5]; [0.5 0.5]; [0.5 0.5]; [0.4 0.6]; ...         % AOCS
                   [0.3 0.7]; [0.3 0.7]; [0.3 0.7]; [0.4 0.6]; [0.4 0.6]; ...                    % TTC
                   [0.8 0.2]; [0.4 0.6]; [0.5 0.5]; [0.3 0.7]; [0.4 0.6] };                      % EPS

problem.maxnfeval = 3e5;

return