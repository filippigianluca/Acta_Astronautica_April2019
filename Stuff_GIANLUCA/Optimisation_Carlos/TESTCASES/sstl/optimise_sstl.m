function [ x, fval, exitflag, output ] = optimise_sstl(problem,par)

    for i = 1:25
        x(i) = 1.0;
    end
    
    for i = 26:31
        x(i) = 0.0;
    end

    fval = problem.objfun(x,problem.par_objfun);
    exitflag = 0;
    output = 0;

return