function [ x, fval, exitflag, output ] = maximise_tof_tweak(problem,par)

    x = nan([1,problem.dim]);
    x(2) = 1.0;
    fval = problem.objfun(x,problem.par_objfun);
    exitflag = 0;
    output = 0;

return