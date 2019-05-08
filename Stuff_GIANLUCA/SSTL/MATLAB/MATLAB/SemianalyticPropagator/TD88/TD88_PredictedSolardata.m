
% %% Predicted data for Kp, Solar F10.7cm radio flux and solar F10.7cm mean radio flux (correct units)
% %   Coveirng the dates of 1/1/2017 to 9/1/2018

load('PredictedDat_2017.mat')

N=2017.00274:0.00259:2018.0219;
%% functions of data for Fx, Fb and kp
solardata.kp=@(d) interp1(N,Kp,d);
solardata.Fx=@(d) interp1(N,F107,d);
solardata.Fb=@(d) interp1(N,F107mean,d);

