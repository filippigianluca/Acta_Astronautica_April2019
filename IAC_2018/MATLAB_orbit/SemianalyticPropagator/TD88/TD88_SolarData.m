
% %% Recored data for Kp, Solar F10.7cm radio flux and solar F10.7cm mean radio flux (correct units)
% %   Coveirng the dates of 1/10/1957 to 1/04/2007

load('Kpindex_data.mat')
load('F107mean_data.mat')
load('F107_data.mat')

N=1957.748:0.00274:2007.28446;
%% functions of data for Fx, Fb and kp
solardata.kp=@(d) interp1(N,Kpindex_data,d);
solardata.Fx=@(d) interp1(N,F107_data,d);
solardata.Fb=@(d) interp1(N,F107mean_data,d);



