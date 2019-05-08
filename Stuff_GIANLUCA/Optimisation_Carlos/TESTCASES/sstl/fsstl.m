function [Costo] = fsstl(d,u,par)

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
    power_needed = (E_day / u(31) + E_night / u(30) )/ daylight;
    eta_1j = u(26)*u(28);
    eta_3j = u(27)*u(29);
    I0 = 1361; % W/m2
    alpha = 123;
    area1 = 0.0007996; %0.0007996; !!!
    area2 = 0.0026639; %0.0007996; !!!p
    C1 = 70/area1;  
    C2 = 300/area2; 

    area_total = power_needed/I0/sin(alpha*pi/180)/(d*eta_1j+(1-d)*eta_3j);

    Costo = area_total*(d*C1+(1-d)*C2);

end
