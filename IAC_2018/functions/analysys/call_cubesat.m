


[d, u] = cubesat_5subsystems_nominal();

interval = 10000:20000;
for i=interval
    d(10) = i; % TIME, fixed (sec)

[P_tot, V] = CUBESAT_5subsystems(d, u, []);

power(i) = P_tot;
data_volume(i) = V;
function_F(i) = V/P_tot;

end

figure
plot(interval, power(interval))

figure
plot(interval, data_volume(interval))

figure
plot(interval, function_F(interval))
