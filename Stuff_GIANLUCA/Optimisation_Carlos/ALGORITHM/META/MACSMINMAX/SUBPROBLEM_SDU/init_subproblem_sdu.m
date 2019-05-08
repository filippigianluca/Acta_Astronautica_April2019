function [subproblem_sdu] = init_subproblem_sdu(problem_minmax)

        % type of problem
        subproblem_sdu.sign_inner = problem_minmax.sign_inner;          % -1 will run minmin

        % objectives
        n_obj = problem_minmax.n_obj;
        subproblem_sdu.n_obj = n_obj;
        for obj = 1:n_obj
            subproblem_sdu.objfun{obj} = @mask_objfun_sdu;            %each function is in the form [f_i] = objfun(i)(d,u,par)
            subproblem_sdu.par_objfun{obj} = struct;                  % the surrogate will get assigned here
            % surrogating
            % These are hardcoded (HC) at the moment but should really come from par_minmax for flexibility
                surrogate.method = 'kriging';
                surrogate.corrfun = @corrgauss;
                surrogate.regrfun = @regpoly0;
                surrogate.training = str2func([lower(surrogate.method) '_training']);
                surrogate.predictor = str2func([lower(surrogate.method) '_predictor']);
                surrogate.model = [];
            subproblem_sdu.par_objfun{obj}.surrogate = surrogate;
        end

        subproblem_sdu.max_set_size = 150;
        subproblem_sdu.set_size_minima = 15; % when building the surrogate it will try to pick this subset of the dataset amongst the best solutions so far.

        % design variables
        dim_d = problem_minmax.dim_d;
        subproblem_sdu.dim_d = dim_d;
        subproblem_sdu.lb_d = zeros(dim_d,1);
        subproblem_sdu.ub_d = ones(dim_d,1);

        % uncertain variables
        dim_u = problem_minmax.dim_u;
        subproblem_sdu.dim_u = dim_u;
        subproblem_sdu.lb_u = cell(1,n_obj);
        subproblem_sdu.ub_u = cell(1,n_obj);
        for obj = 1:n_obj
            subproblem_sdu.lb_u{obj} = repmat({0},[dim_u,1]);
            subproblem_sdu.ub_u{obj} = repmat({1},[dim_u,1]);
        end

        % % maxnfeval: optional to give it here or in the algorithm
        % problem_minmax.maxnfeval = 5.0e3;
        % problem_minmax.nfouter_nested = 50;

        

return
