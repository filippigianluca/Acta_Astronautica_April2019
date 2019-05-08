function [masked] = mask_objfun_minmarek_outer(d,par_objfun)

[y0,mse0] = par_objfun.surrogate.predictor(d, par_objfun.surrogate.model);
y = y0(end-par_objfun.n_obj+1:end);
mse = mse0(end-par_objfun.n_obj+1:end);

masked = -par_objfun.surrogate.indicator(y,mse,par_objfun.ymin); % negative because we will always maximize PI or EI

return