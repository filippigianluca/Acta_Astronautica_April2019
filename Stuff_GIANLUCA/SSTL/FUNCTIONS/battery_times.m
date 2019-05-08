function t_on_off = battery_times(orb_params,MJD)

% orb_ params = [a e i RAAN omega theta tf dt] for every LAE burn.
% the 1st row is injection by launcher. the parameters are after the burn. tf and dt in hours after launch
	mu = astro_constants(13);
	R = astro_constants(23);


	t_on_off = [];
	n_burns = size(orb_params,1);
	t_after_last_burn = 24;
	in_eclipse_start=false;
	in_eclipse_end = false;
	for b=1:n_burns
		params_b = orb_params(b,1:6);
		tf_b = orb_params(b,7);
		tt = MJD+tf_b/24.0;
		dt_b = orb_params(b,8);
		
		if b == n_burns
			tf_arc = tf_b + t_after_last_burn;
		else
			tf_arc = orb_params(b+1,7)-orb_params(b+1,8);
		end

		Eq = kep2eq(params_b);
		[L_i, L_o, kep, dt, kep2, alpha_Sun, beta_Sun, phi_Sun] = sol_Eclipse_Eq(Eq, tt, mu, R);
      
		L = Eq(6);
		if ~isnan(L_i)
			dLi = mod(L-L_i, 2*pi);
			dLe = mod(L_o-L_i, 2*pi);
			in_eclipse_start = 0<= dLi && dLi <= dLe;
			if in_eclipse_start
% 				dt_off = kepEq_t(L_o -kep(4) - kep(5), kep2(1), kep2(2), mu, L - kep(4) - kep(5), 0);
                dt_off = kepEq_t_matlab(L_o -kep(4) - kep(5), kep2(1), kep2(2), mu, L - kep(4) - kep(5), 0);
				dt_off = dt_off +  2*pi*sqrt(kep(1)^3/mu)*(dt_off <0);
				dt_off = dt_off/3600.0;
				if ~in_eclipse_end % otherwise, there is a last t_on but no last t_off
					t_on_off(end+1,1) = tf_b - dt_b; % even if dt_b == 0 because we are in eclipse
				end
				t_off = tf_b + dt_off;
				if t_off < tf_arc || b==n_burns
					t_on_off(end,2) = t_off;
					in_eclipse_end = false;
				else
					in_eclipse_end = true;
				end
			else
				% first account for the LAE firing "eclipse"
				if (in_eclipse_end)
					t_on_off(end,2) = tf_b;
				elseif dt_b > 0
					t_on_off(end+1,1) = tf_b -dt_b;
					t_on_off(end,2) = tf_b;
				end

				% then account for a potential next eclipse
% 				dt_on = kepEq_t(L_i -kep(4) - kep(5), kep2(1), kep2(2), mu, L - kep(4) - kep(5), 0);
                dt_on = kepEq_t_matlab(L_i -kep(4) - kep(5), kep2(1), kep2(2), mu, L - kep(4) - kep(5), 0);
				dt_on = dt_on +  2*pi*sqrt(kep(1)^3/mu)*(dt_on <0);
				dt_on = dt_on/3600.0;
				t_on = tf_b + dt_on;
				t_off = t_on + dt/3600.0;
				if t_on < tf_arc %otherwise nothing to do
					t_on_off(end+1,1) = t_on;
					if t_off < tf_arc || b == n_burns
						t_on_off(end,2) = t_off;
						in_eclipse_end = false;
					else
						in_eclipse_end = true;
                    end
				end
			end
		else % there is no eclipse in the current orbit
			if (in_eclipse_end)
				t_on_off(end,2) = tf_b;
				in_eclipse_end = false;
			elseif dt_b > 0
				t_on_off(end+1,1) = tf_b -dt_b;
				t_on_off(end,2) = tf_b;
			end
		end
	end
	if (in_eclipse_end)
		t_on_off(end,2) = t_off;
	end

end
