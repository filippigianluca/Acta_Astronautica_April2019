function [minmax] = evaluate_minmax(problem_minmax, algo_minmax, algo_outer, algo_inner)

%% evaluate the worst-case scenario (minmax)



global nfevalglobal;
nfevalglobal = 0;



problem_minmax.sign_inner = 1;  %  for minmax


%% optimise
[ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax); % algo_minmax.optimise



minmax.d = dmin;
minmax.u = output.u;
minmax.f = fminmax;

end