function J = cost_fun(ctrl, m0, Isp, g0, L_star, T_star)

% NOTE: delta_v is non-dimensional


% Tsiolkovsky equation
% delta_v = Isp*g0*log(m0/mf);

% delta_m = m0*(1-exp(-delta_v/Isp/g0))
J = (m0*(1 -exp(-sum(ctrl(2:4:end)*T_star/L_star +ctrl(3:4:end)*T_star/L_star +ctrl(4:4:end)*T_star/L_star) /Isp/g0)))^2;

