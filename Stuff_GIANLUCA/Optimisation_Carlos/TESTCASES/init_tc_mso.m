function [problem_minmax] = init_tc(varargin)

if (nargin == 0)
    id = 1;
else
    id = varargin{1};
end


switch (id)
    case 11
        %% initialise problem TC1_F1 (MV1)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
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
        problem_minmax.maxnfeval = 1e4*dim_d; % 1.0e3;

    case 12
        %% initialise problem TC1_F2 (MV3)

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
        problem_minmax.maxnfeval = 1e4*dim_d; % 1.0e3;
        
    case 21

        %% initialise problem TC2_F1 (MV2)

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
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e4;

    case 22

        %% initialise problem TC2_F2 (MV8)

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
        problem_minmax.maxnfeval = 1e4*dim_d; % 1e4;


    case 31
        %% initialise problem TC3_F1 (MV2)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv2};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e4;
        
    case 32
        %% initialise problem TC3_F2 (EM1)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 1.0e4;

    case 41
        %% initialise problem TC4_F1 (MV8)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 1e4;

    case 42
        %% initialise problem TC4_F2 (MV9)

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
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e3;

    case 51
        %% initialise problem TC5_F1 (MV8)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 1e4;

    case 52
        %% initialise problem TC5_F2 (EM1)

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
        problem_minmax.maxnfeval = 1e4*dim_d; % 1.0e4;

    case 61
        %% initialise problem TC6_F1 (MV10)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv10};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 1;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-4,[dim_d,1]);
        problem_minmax.ub_d = repmat(2*pi,[dim_d,1]);

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[pi,5,5.5]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[4,6,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e3;
    
    case 62
        %% initialise problem TC6_F2 (MV9)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv9};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 1;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-4,[dim_d,1]);
        problem_minmax.ub_d = repmat(2*pi,[dim_d,1]);

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-pi/2,0,3*pi/4]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-pi/6,pi,3*pi/2]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e3;
    
    case 71

        %% initialise problem TC7_F1 (MV2)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv2};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 5.0e4;

    case 72

        %% initialise problem TC7_F2 (MV8)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 1e4;

    case 73
        %% initialise problem TC7_F3 (EM1)

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1e4*dim_d; % 1.0e4;

    case 20
        %% initialise problem TC1_scalarised

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 1;
        problem_minmax.objfun = {@mv_1_3};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct};
        problem_minmax.par_objfun{1}.lambda = [.5 .5];

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
        problem_minmax.maxnfeval = 2e4*dim_d; % 1.0e3;

    otherwise
        error('init_testcase: trying to initialise testcase with unknown id')
end
return