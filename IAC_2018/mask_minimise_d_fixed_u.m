function varargout = mask_minimise_d_fixed_u(d,par)



u = par{1, 1}.fix.nominal_u;
F = par{1, 1}.objfun{1}(d, u, par{1});


% objfun and coinstraint can be in the same file
if isstruct(F)
    varargout{1} = F.f;
    varargout{2} = F.c;
    varargout{3} = F.ceq;
else
    varargout{1} = F;
end
end