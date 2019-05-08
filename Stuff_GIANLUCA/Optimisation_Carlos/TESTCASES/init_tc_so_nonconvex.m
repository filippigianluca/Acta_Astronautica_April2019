function [problem_minmax] = init_tc_so_convex(varargin)

if (nargin == 0)
    id = 1;
else
    id = varargin{1};
end


switch (id)
    case 1
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 16;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 16;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 5.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 2e3;

    case 2
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv2};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 5.0e3;
        % problem_minmax.nfouter_nested = 100;
        problem_minmax.maxnfeval_eff = 1e3;

    case 3
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv3};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 15.0e3;
        % problem_minmax.nfouter_nested = 200;
        problem_minmax.maxnfeval_eff = 4e3;
        
    case 4
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 7.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 1e5;

    case 5
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv9};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-pi/2,0,3*pi/4]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-pi/6,pi,3*pi/2]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 4.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff =2e4;

    case 6
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv10};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-4.0,[dim_d,1]);
        problem_minmax.ub_d = repmat(2*pi,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[pi,5,5.5]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[4,6,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 50.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 3e4;

    case 7
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 90.0e3;
        % problem_minmax.nfouter_nested = 100;
        problem_minmax.maxnfeval_eff = 8e3;

    case 8
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp10};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 1;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1e-4,[dim_d,1]);
        problem_minmax.ub_d = repmat(10,[dim_d,1]);

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({1e-4},[dim_u,1])};
        problem_minmax.ub_u = {repmat({10},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 5.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 1e3;

    case 9
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp11};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 1;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(10,[dim_d,1]);

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({10},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 5.0e3;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 2e3;
        
    case 10
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_rbf_2d};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 2.0e3;
        % problem_minmax.nfouter_nested = 100;
        problem_minmax.maxnfeval_eff = 2e3;
    case 11
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_rbf_5d};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 5;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(10,[dim_d,1]);

        % uncertain variables
        dim_u = 5;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 10.0e3;
        % problem_minmax.nfouter_nested = 200;
        problem_minmax.maxnfeval_eff = 6e3;

    case 12
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_rbf_10d};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 10;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(10,[dim_d,1]);

        % uncertain variables
        dim_u = 10;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 1.0e5;
        % problem_minmax.nfouter_nested = 250;
        problem_minmax.maxnfeval_eff = 5e4;

    case 13
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_sonc_13};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 6;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-6,[dim_d,1]);
        problem_minmax.ub_d = repmat(6,[dim_d,1]);

        % uncertain variables
        dim_u = 6;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-5.12},[dim_u,1])};
        problem_minmax.ub_u = {repmat({5.12},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 1.0e5;
        % problem_minmax.nfouter_nested = 250;
        problem_minmax.maxnfeval_eff = 2e4;

    case 14
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_sonc_14};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 25;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(1,[dim_d,1]);

        % uncertain variables
        dim_u = 25;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 1.0e5;
        % problem_minmax.nfouter_nested = 50;
        problem_minmax.maxnfeval_eff = 1e4;
    case 15
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_sonc_15};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 10;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 10;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 1.0e5;
        % problem_minmax.nfouter_nested = 250;
        problem_minmax.maxnfeval_eff = 1e5;

    case 16
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_damper};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = [0 ; 0];
        problem_minmax.ub_d = [1 ; 2];

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0.0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({2.5},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 1.0e5;
        % problem_minmax.nfouter_nested = 250;
        problem_minmax.maxnfeval_eff = 8e2;
    otherwise
        error('init_testcase_so: trying to initialise testcase with unknown id')
end
return