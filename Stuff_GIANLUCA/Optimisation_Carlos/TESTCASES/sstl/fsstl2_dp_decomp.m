function [dp] = fsstl2_dp_decomp(d,u,par)

    daylight = 63 * 60 * length(par.day);
    computed_orbits = nan(length(par.orbits),2);
    E_day = 0;
    E_night = 0;

    for id_orbit = par.day
        if isnan(computed_orbits(id_orbit,1))
            orbit = par.orbits{id_orbit};
            computed_orbits(id_orbit,:) = [0 0];
            orbit_time = 0;
            for i = 1:size(orbit,1)
                if orbit(i,2) > 0
                    timestep = orbit(i,2);
                    next_orbit_time = orbit_time + timestep;
                    if (next_orbit_time <= 19*60 || orbit_time > (19+63)*60) % all night
                        timestep_night = timestep;
                        timestep_day = 0;
                    elseif (orbit_time > 19*60 && next_orbit_time <= (19+63)*60) % all daylight
                        timestep_day = timestep;
                        timestep_night = 0;
                    elseif (orbit_time <= 19*60) % begins night ends day or night
                        timestep_night = 19*60-orbit_time;
                        timestep_day = next_orbit_time - 19*60;
                        timestep_night = max (timestep_night, timestep_night + timestep_day - 63*60);
                        timestep_day = min(timestep_day, 63*60);
                        % if (timestep_day > 63 * 60)
                        %     timestep_night = timestep_night + timestep_day - 63*60;
                        %     timestep_day = 63*60;
                        % end
                    else % begins day ends night
                        timestep_day = (19+63)*60 - orbit_time;
                        timestep_night = next_orbit_time - (19+63)*60;
                    end
                    power = 0;
                    for j = 1:20
                        if (orbit(i,2+j)>0)
                            powers_to_sum = par.groups{j};
                            power = power + orbit(i,2+j)*sum(u(powers_to_sum))/100;
                        end
                    end
                    % power
                    % timestep_day
                    % timestep_night
                    % a=[]
                    orbit_time = next_orbit_time;
                    computed_orbits(id_orbit,:) = computed_orbits(id_orbit,:) + [power * timestep_day, power * timestep_night];
                end
            end
        end

        E_day = E_day + computed_orbits(id_orbit,1);
        E_night = E_night + computed_orbits(id_orbit,2);

    end

    % real_power = (E_day + E_night)/82726
    power_needed = (E_day / (1.0-u(31)) + E_night / (1.0-u(30)) )/ daylight;
    eta_1j = (1.0-u(26))*(1.0-u(28));
    eta_3j = (1.0-u(27))*(1.0-u(29));

    I0 = 1361; % W/m2
    alpha = 123;
    mu = d(1);
    A = d(2);

    dp = power_needed - I0*sin(alpha*pi/180)*A*(mu*eta_1j+(1-mu)*eta_3j);

end
