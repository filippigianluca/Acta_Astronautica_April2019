function [problem_minmax] = init_problem_transfer()

    %% initialise problem TC1

    % type of problem
    problem_minmax.sign_inner = 1;          % -1 will run minmin

    % objectives
    problem_minmax.n_obj = 2;
    problem_minmax.objfun = {@objfun_tof,@objfun_dv_transfer};   %each function is in the form [f_i] = objfun(i)(d,u,par)
    problem_minmax.par_objfun = {struct,struct};

    % parameters for dv function
    par = struct;
    % ephemeris [a e i W w M] in AU and degrees
    % %Ceres
    % ephref0 = [2.767046248500289    .07553461024389638  10.5935097971363    80.30991865594387   73.11534200131032    352.2304611765882];
    % sigma0  = [2.6002e-11           1.1644e-11          1.9258e-09          1.1728e-08          1.1437e-08           2.6858e-09];    
    % % Pallas
    % ephref1 = [2.772766124346009    .2305054995905846   34.83687567654144   173.0838008845752   310.0063356539786    334.3231878791698];
    % sigma1  = [4.3096e-09           3.3819e-08          3.7277e-06          6.4667e-06          9.5268e-06           9.1161e-06];
    % Pirogov
    ephref0 = [2.89900113535035     .01886105155538726  2.165374950008691   164.874796730799    283.3147692190075   54.13276106257043];
    sigma0  = [1.4206e-08           3.4064e-08          4.0239e-06          9.9414e-05          0.00014519          0.00010605];
    % Souten
    ephref1 = [2.622105895759311    .1259482484043224    1.605515506620147  156.9586376277082   296.4015148476117    269.1495992707831];
    sigma1  = [1.1948e-08           3.4848e-08           4.2811e-06         0.00014775          0.00014887           1.9444e-05];

    % ephemeris adim
    ephref0(3:6) = deg2rad(ephref0(3:6));
    ephref1(3:6) = deg2rad(ephref1(3:6));
    sigma0(3:6) = deg2rad(sigma0(3:6));
    sigma1(3:6) = deg2rad(sigma1(3:6));
    par.ephref0         = ephref0;        
    par.ephref1         = ephref1;    
    par.T0              = 2*pi*sqrt(ephref0(1)^3);;
    par.T1              = 2*pi*sqrt(ephref1(1)^3);
    par.epoch_ref       = 58200;

    mu = 1.32712440018e11;
    au = 1.49597870691e8;
    tu = sqrt(au^3/mu);
    tu2day = tu/86400;
    day2tu = 1/tu2day;

    par.au              = au;
    par.tu              = tu;
    par.day2tu          = day2tu;    

    problem_minmax.par_objfun{2} = par;
    % Design variables
    dim_d = 2;
    problem_minmax.dim_d = dim_d;
    problem_minmax.lb_d = [61500; 1.5];
    problem_minmax.ub_d = [63500; 500];

    % Uncertain variables
    dim_u = 14;
    problem_minmax.dim_u = dim_u;
    k_sigma = 4;
    % bounds
    lb_u = cell(dim_u,1);
    ub_u = cell(dim_u,1);
    bpa  = cell(dim_u,1);
    lb_u{1,1} = [-1 0 1 4];     ub_u{1,1} = [0 1 4 7];      bpa{1,1}  = [0.25 0.25 0.25 0.25];
    lb_u{2,1} = [-1 0 1 3];     ub_u{2,1} = [0 1 3 5];      bpa{2,1}  = [0.25 0.25 0.25 0.25];

    for i=3:8
        lb_u{i,1} = [-k_sigma*sigma0(i-2) 0];
        ub_u{i,1} = [0 +k_sigma*sigma0(i-2)];
        bpa{i,1}  = [0.5 0.5];
    end
    for i=9:14
        lb_u{i,1} = [-k_sigma*sigma1(i-8) 0];
        ub_u{i,1} = [0 +k_sigma*sigma1(i-8)];
        bpa{i,1}  = [0.5 0.5];
    end
    problem_minmax.lb_u = {lb_u,lb_u};
    problem_minmax.ub_u = {ub_u,ub_u};
    problem_minmax.bpa  = {bpa , bpa};
    % maxnfeval: optional to give it here or in the algorithm
    problem_minmax.maxnfeval = 5e5;

 
return