function [record,d_info,u_record,nfeval_cc] = cross_check(problem_fix_d,record,d_info,u_record,cc_opt)

	sign_inner = problem_fix_d.par_objfun.sign;
	nfeval_cc = 0;

	if (cc_opt.to_crosscheck > 0) % otherwise nothing to do
		u_record_aux = cell(1,cc_opt.n_obj);
		func = problem_fix_d.objfun;

		n_d = cc_opt.dim_d;
		n_u = cc_opt.dim_u;
		n_obj = cc_opt.n_obj;
		skip = [];
		if isfield(cc_opt,'skip')
			skip = cc_opt.objectives_to_skip
		end

		i_to_check = [];
		i_checked = [];
		if (cc_opt.to_crosscheck > 1)
			i_to_check = 1:size(record,1);
		else
			i_to_check = find(dominance(record(:,end-n_obj+1:end),0) == 0);
			if(n_obj>1) 
            	i_to_check = i_to_check';
        	end
		end

		while(~isempty(i_to_check))
			% cross-check who needs it
			for i=i_to_check
				problem_fix_d.par_objfun.d = record(i,1:n_d);
				for obj=1:n_obj
					if ~ismember(obj,skip)
						problem_fix_d.par_objfun.objective = obj;
						for j = 1:size(u_record{obj},1)
							u_du = u_record{obj}(j,:);
							if ~ismember(j,d_info(i).u_checked{obj})
				                [f_du] = func (u_du, problem_fix_d.par_objfun);
				                nfeval_cc = nfeval_cc + 1;
				                d_info(i).u_checked{obj}(end+1) = j;
				                if f_du < -sign_inner*record(i,n_d+n_u*n_obj+obj)
				                	record(i,n_d+n_u*n_obj+obj) = -sign_inner*f_du;
				                	record(i,n_d+(obj-1)*n_u+1 : n_d+obj*n_u) = u_du;
				                	d_info(i).umax(obj) = j;
				                	d_info(i).fmax(obj) = -sign_inner*f_du;
				                end
				            end
						end
					end
				end
			end

			if cc_opt.ls_flag_crosscheck > 0 % you finished cross checking so you do the local searches
				for i=i_to_check
					problem_fix_d.par_objfun.d = record(i,1:n_d);
					for obj = 1:n_obj
						if ~ismember(obj,skip)
							problem_fix_d.par_objfun.objective = obj;
							to_ls = [];
							if cc_opt.ls_flag_crosscheck == 1
								to_ls =  d_info(i).umax(obj); % only from the max of the cross-checks
							else
								to_ls = 1:size(u_record{obj},1); % from everyone in u_record
							end

							l_aux = size(u_record_aux{obj},1)+1; % if u_ls maximises, it will be allocated here in u_record_aux
							l = l_aux + size(u_record{obj},1);	 % and this will be its position in u_record when it is updated

							for j = to_ls
								if ~ismember(j,d_info(i).u_ls{obj}) % run local search unless you already did. otherwise nothing to do
									cc_opt.algo_ls{obj}.par.initial_population = u_record{obj}(j,:);
									[ u_ls, f_ls , ~ , output_aux] = cc_opt.algo_ls{obj}.optimise(problem_fix_d,cc_opt.algo_ls{obj}.par);
									nfeval_cc = nfeval_cc + output_aux.funcCount;
									d_info(i).u_ls{obj}(end+1) = j;

									[~, loc1] = ismember(round(cc_opt.tol*u_ls),round(cc_opt.tol*u_record{obj}),'rows');
									if loc1 > 0 % you already have this u in the archive. no need to add it
					                    % but maybe you had never checked this d and this u together
					                    if ~ismember(loc1,d_info(i).u_checked{obj})
					                        d_info(i).u_checked{obj}(end+1) = loc1;
					                    end
					                    % NOTE: I ommit the analogous check for ls because I think the output of an ls
					                    % is better not marked in u_ls in case the search didn't converge.

					                    % and maybe it maximises
					                    if f_ls < -sign_inner*record(i,n_d+n_u*n_obj+obj)
					                        record(i , n_d+(obj-1)*n_u+1 :n_d+ obj*n_u) = u_ls;
					                        record(i , n_d+n_u*n_obj+obj) = -sign_inner*f_ls;
					                        d_info(i).umax(obj) = loc1;
					                        d_info(i).fmax(obj) = -sign_inner*f_ls;
					                    end

					                else % there is a new u associated that is not in u_record
					                	% but is it in u_record_aux?
	                                    loc1_aux = 0;
	                                    if ~isempty(u_record_aux{obj})
	                                        [~, loc1_aux] = ismember(round(cc_opt.tol*u_ls),round(cc_opt.tol*u_record_aux{obj}),'rows');
	                                    end
					                	if loc1_aux > 0
					           				% it is in u_record_aux, so no need to add it
					                		% but maybe you had never checked this d and this u together (for sure)
					                		true_loc1_aux = loc1_aux + size(u_record{obj},1);
						                    if ~ismember(true_loc1_aux,d_info(i).u_checked{obj})
						                        d_info(i).u_checked{obj}(end+1) = true_loc1_aux;
						                    end
						                    % NOTE: I ommit the analogous check for ls because I think the output of an ls
						                    % is better not marked in u_ls in case the search didn't converge.

						                    % and maybe it maximises
						                    if f_ls < -sign_inner*record(i,n_d+n_u*n_obj+obj)
						                        record(i , n_d+(obj-1)*n_u+1 :n_d+ obj*n_u) = u_ls;
						                        record(i , n_d+n_u*n_obj+obj) = -sign_inner*f_ls;
						                        d_info(i).umax(obj) = true_loc1_aux;
						                        d_info(i).fmax(obj) = -sign_inner*f_ls;
						                    end


					                	else %it is completely new

					                		%% actually no need to store local minima. so this below is commented out

					                		% u_record_aux{obj}(end+1,:) = u_ls;     % archive it
						                    % l = size(u_record{obj},1)+size(u_record_aux{obj},1);       
						                    % d_info(i).u_checked{obj}(end+1) = l;   % remember it has already been cross-checked
					                		

					                		% now does it maximise?
						                    if f_ls < -sign_inner*record(i,n_d+n_u*n_obj+obj)
						                    	% then archive it
						                    	u_record_aux{obj}(l_aux,:) = u_ls;     % archive it overwriting local minima
						                		if ~ismember(l,d_info(i).u_checked{obj})
							                    	d_info(i).u_checked{obj}(end+1) = l;   % remember it has already been cross-checked
						                		end
						                		% and do the usual stuff
						                        record(i , n_d+(obj-1)*n_u+1 :n_d+ obj*n_u) = u_ls;
						                        record(i , n_d+n_u*n_obj+obj) = -sign_inner*f_ls;
						                        d_info(i).umax(obj) = l;
						                        d_info(i).fmax(obj) = -sign_inner*f_ls;

						                    end     
					                	end
					                end
								end
							end
						end
					end
				end
			end

			i_checked = [i_checked, i_to_check];

			if (cc_opt.to_crosscheck > 1)
				i_to_check = [];
			else
				i_to_check = find(dominance(record(:,end-n_obj+1:end),0) == 0);
				if(n_obj>1) 
            		i_to_check = i_to_check';
        		end
        		i_to_check = i_to_check(~ismember(i_to_check,i_checked));
            end
            
            
		end


		% append the u_record_aux to u_record
		for obj = 1:n_obj
			if ~ismember(obj,skip)
				u_record{obj} = [u_record{obj}; u_record_aux{obj}];
			end
		end
	end




return