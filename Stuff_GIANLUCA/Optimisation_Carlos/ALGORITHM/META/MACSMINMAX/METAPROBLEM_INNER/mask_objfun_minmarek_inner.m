function [masked] = mask_objfun_minmarek_inner(u,par_objfun)

[y,mse] = par_objfun.surrogate.predictor(u, par_objfun.surrogate.model);
if par_objfun.use_indicator
	masked = -par_objfun.surrogate.indicator(y,mse,par_objfun.ymin); % negative because we will always maximize PI or EI
else
	masked = y;
end

return