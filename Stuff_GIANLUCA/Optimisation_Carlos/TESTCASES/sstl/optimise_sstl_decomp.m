function [ x, fval, exitflag, output ] = optimise_sstl_decomp(problem,par)

    x = problem.ub;

    fval = problem.objfun(x,problem.par_objfun);
    exitflag = 0;
    output = 0;

return