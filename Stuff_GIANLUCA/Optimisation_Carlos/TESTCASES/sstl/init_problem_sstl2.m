function [problem] = init_problem_sstl2()

% type of problem
problem.sign_inner = 1; % -1 will run minmin
% objectives
problem.n_obj = 1;
problem.objfun = {@fsstl2_dp};

%% DESIGN VECTOR d

problem.dim_d = 2;
problem.lb_d = [0 ; 0];
problem.ub_d = [1 ; 0.36];

% different exchange
problem.dim_u = 31;
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
    problem.lb_u{1}{i,1} = [info_u(i,1), info_u(i,1)*(1+0.5*info_u(i,2))];
    problem.ub_u{1}{i,1} = [info_u(i,1)*(1+0.5*info_u(i,2)), info_u(i,1)*(1+info_u(i,2))];

    if info_u(i,2) == 0.02
        problem.lb_u{1}{i,1} = [info_u(i,1)];
        problem.ub_u{1}{i,1} = [info_u(i,1)*(1+info_u(i,2))];    
    end

end
problem.lb_u{1}{i+1,1} = [16.42,19]/100; % efficiency single
problem.ub_u{1}{i+1,1} = [19,21.16]/100;

problem.lb_u{1}{i+2,1} = [21.7,25]/100; % efficiency triple
problem.ub_u{1}{i+2,1} = [25,29.81]/100;

problem.lb_u{1}{i+3,1} = [0.92,0.95]; % assembly single
problem.ub_u{1}{i+3,1} = [0.95,0.99];

problem.lb_u{1}{i+4,1} = [0.93,0.95]; % assembly triple
problem.ub_u{1}{i+4,1} = [0.95,0.97];

problem.lb_u{1}{i+5,1} = [0.94,0.96]; % efficiency night
problem.ub_u{1}{i+5,1} = [0.96,0.98];

problem.lb_u{1}{i+6,1} = [0.95,0.97]; % efficiency day
problem.ub_u{1}{i+6,1} = [0.97,0.99];

problem.bpa{1} = cell(size(problem.lb_u{1},1),1);
for i = 1:size(problem.lb_u{1},1)
    n_int = size(problem.lb_u{1}{i,1},2);
    problem.bpa{1}{i,1} = repmat(1/n_int,[1, n_int]);
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


return