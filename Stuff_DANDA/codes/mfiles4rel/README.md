# functions for reliability calculation
interface through reliability_SYSTEM.m

parameters (of reliabnility_SYSTEM):
- t, time for which we want to calculate reliability (i.e. probability that system will function from t=0 to t=t)
- d, vector of design parameters, all the design parameters for the subsystems are stacked here, i.e. first n_aocs will be design parameters for AOCS, next n_power will be design variables of power, etc...
- u, vector of uncertain variables, same stacking as with design vector
