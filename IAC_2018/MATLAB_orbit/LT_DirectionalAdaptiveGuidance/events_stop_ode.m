function [value,isterminal,direction] = events_stop_ode(t, x, DAG_parameters)


error = (DAG_parameters.target_flag' .* abs(DAG_parameters.x_target(1:5)' - x(1:5)));

value = 1000;

if error(1) < 1e-3 && error(2) < 1e-2 && error(3) < 1e-3
    value = 0;
end

isterminal = 1;

direction = 0;



return





