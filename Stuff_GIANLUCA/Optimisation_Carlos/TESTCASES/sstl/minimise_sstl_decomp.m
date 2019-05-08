function [ x, fval, exitflag, output ] = minimise_sstl_decomp(problem,par)

    x = problem.lb;

    fval = problem.objfun(x,problem.par_objfun);
    exitflag = 0;
    output = 0;

return