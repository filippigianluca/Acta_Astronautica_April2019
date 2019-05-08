% Optimisation of low-thrust transfer.
% Direct optimisation with multiple shooting and alternation of coast and
% thrust arc. Analytic propagation of the coast arcs.

%% LT_transfer_main.m

% Low-thrust transfer given initial state, final state and time of flight.
% It is possible to solve the problem for feasibility of to minimise the
% DeltaV

%% LT_transfer_main2.m

% Low-thrust transfer given states of departure and arrival body at the
% same initial epoch (initial time). 
% It is possible to solve the problem:
% - feasibility
% - optimisation of time of fligth
% - optimisation of DeltaV
% The difference with LT_transfer_main.m is that in this case the time of
% flight is not given 