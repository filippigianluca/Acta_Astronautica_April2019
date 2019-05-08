%---------------------------------------------------------------------------
%
% TrjData:
%   The TrjDate class provides a framework for the integration of spacecraft
%   trajectories including both the state and the variational equations. 
%   It handles the required initialization steps, the packing and unpacking
%   of the state vector, state transition matrix and sensitivity matrix and
%   the mapping between the external time scale (Modified Julian Date TT)
%   and the internal time scale (seconds since reference epoch)
%
% Last modified:   2015/08/12   M. Mahooti
%
%---------------------------------------------------------------------------
function [y,Y] = Trj(n_var,Mjd_UTC,y,Y)

global AuxParam

% Initialization
t = (Mjd_UTC-AuxParam.Mjd_UTC)*86400;

% Internal integration to given time
options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
[~,yout] = radau(@Accel,[0 t],y,options);
y = yout(end,:)';
Y = DEInteg(@VarEqn_,0,t,1e-13,1e-6,7*n_var,Y);

