function [f] = mask_obfjun_sdu(d,u,par)
    [y,mse] = par.surrogate.predictor([d, u], par.surrogate.model);
    f=y;
return