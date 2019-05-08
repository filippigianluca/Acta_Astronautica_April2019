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
        problem_minmax.objfun = {@mwp1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-5},[dim_u,1])};
        problem_minmax.ub_u = {repmat({5},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 5.0e3;
        problem_minmax.nfouter_nested = 50;

    case 2
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp2};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-5},[dim_u,1])};
        problem_minmax.ub_u = {repmat({5},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 5.0e3;
        problem_minmax.nfouter_nested = 100;

    case 3
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp3};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-3},[dim_u,1])};
        problem_minmax.ub_u = {repmat({3},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 15.0e3;
        problem_minmax.nfouter_nested = 200;
        
    case 4
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp4};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 3;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-3},[dim_u,1])};
        problem_minmax.ub_u = {repmat({3},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 7.0e3;
        problem_minmax.nfouter_nested = 50;
        
    case 5
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp5};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 3;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 3;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 4.0e3;
        problem_minmax.nfouter_nested = 50;
        
    case 6
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp6};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 3;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-2},[dim_u,1])};
        problem_minmax.ub_u = {repmat({2},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 50.0e3;
        problem_minmax.nfouter_nested = 50;
        
    case 7
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp7};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 5;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 5;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-3},[dim_u,1])};
        problem_minmax.ub_u = {repmat({3},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 90.0e3;
        problem_minmax.nfouter_nested = 100;
        
    case 8
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
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
        problem_minmax.maxnfeval = 0.5e3;
        problem_minmax.nfouter_nested = 100;
        
    case 9
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp9};            %each function is in the form [f_i] = objfun(i)(d,u,par)
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
        problem_minmax.maxnfeval = 4.0e3;
        problem_minmax.nfouter_nested = 60;
        
        
    case 10
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp12};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = [-0.5 ; 0 ];
        problem_minmax.ub_d = [ 0.5 ; 1 ];

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({10},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 2.0e3;
        problem_minmax.nfouter_nested = 100;
        
    case 11
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mwp13};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-1,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({10},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 10.0e3;
        problem_minmax.nfouter_nested = 200;
            
    case 12
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_soc_12};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(25,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({0},[dim_u,1])};
        problem_minmax.ub_u = {repmat({25},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e5;
        problem_minmax.nfouter_nested = 1.25e3;
    
    case 13
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_soc_13};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};
        
        % design variables
        dim_d = 6;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5.12,[dim_d,1]);
        problem_minmax.ub_d = repmat(5.12,[dim_d,1]);

        % uncertain variables
        dim_u = 6;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-6},[dim_u,1])};
        problem_minmax.ub_u = {repmat({6},[dim_u,1])};


        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e5;
        problem_minmax.nfouter_nested = 1e4;

    case 14
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_soc_14};            %each function is in the form [f_i] = objfun(i)(d,u,par)
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
        problem_minmax.maxnfeval = 1.0e5;
        problem_minmax.nfouter_nested = 50;
    
    case 15
        %% initialise problem

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@f_soc_15};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 16;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-5,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 16;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({-1},[dim_u,1])};
        problem_minmax.ub_u = {repmat({1},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e5;
        problem_minmax.nfouter_nested = 1.5e3;
    otherwise
        error('init_testcase_so: trying to initialise testcase with unknown id')
end
return