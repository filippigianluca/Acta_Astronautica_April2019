function [masked] = mask_objfun_mo_sa_outer(d,par_objfun)

[y,~] = par_objfun.surrogate.predictor(d, par_objfun.surrogate.model);

masked = y(end-par_objfun.n_obj+1:end);

return