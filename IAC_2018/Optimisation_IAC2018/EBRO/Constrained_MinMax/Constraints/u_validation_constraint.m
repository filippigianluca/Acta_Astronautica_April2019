% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [ARCHIVE_crosscheck, nfeval,all_f, violation] = u_validation_constraint(problem_fix_d, c_archive, local_search_flag,objectives, ARCHIVE_normalised, ARCHIVE_crosscheck)
        
violation   = [];
violation_C = [];

minmaxsign = problem_fix_d.par_objfun.sign;
archsize = size(ARCHIVE_crosscheck{1, 1},1);
u = cell(1,max(objectives));
f = nan(archsize,max(objectives));
for obj = objectives
    u{obj} = nan(archsize,problem_fix_d.dim);
end
nfeval = 0;
func = problem_fix_d.objfun;
all_f=[];
if nargout > 3
    all_f = cell(1,max(objectives)); % all evaluations stored like all_f{obj}(idx_d_record,idx_u_record)
    % careful when local search, we are not storing u_f only u_0 and it's f(d,u_f)
end

for idx_d = 2:archsize
    % fix d and scale
    problem_fix_d.par_objfun.d = ARCHIVE_normalised{1, 1}{idx_d, 1};  %d_record(idx_d,:);
    for obj = objectives
        
        map_info = problem_fix_d.par_objfun.map_u_info{obj};
        problem_fix_d.par_objfun.objective = obj;
        f_du = [];
        u_du=[];
        f_d= - ARCHIVE_normalised{1, 1}{idx_d, 4}; %-f_archive(idx_d);
        u_d= ARCHIVE_normalised{1, 1}{idx_d, 2}; %nan(1,problem_fix_d.dim);
        
        for idx_u = 2:size(ARCHIVE_normalised{1, 1},1) %1:size(u_record{obj},1)
            u_0=ARCHIVE_normalised{1, 1}{idx_u, 2};  %u_record{obj}(idx_u,:);
            c_du = [];
            ceq_du = [];
            if (local_search_flag)
                stop = 0;
                iter = 0;
                crowding = 0.1;
                lb_local = u_0 - crowding/2;
                ub_local = u_0 + crowding/2;
                lb_local(lb_local < 0) = 0;
                ub_local(ub_local > 1) = 1;
                options = optimset('Display','none','MaxFunEvals',50*problem_fix_d.dim,'TolFun',1e-8,...%'LargeScale','off',...
                    'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
                
                %----------------------------------------------------------
                if ~isempty(problem_fix_d.fitnessfcn.constr)                    
                    % CONSTRAINT
%                     [constraint] = feval(problem_fix_d.par_objfun.mask_constraints, u_du, problem_fix_d.par_objfun);
                    func_constr = problem_fix_d.par_objfun.mask_constraint;
                else
                    func_constr = [];
                end
                %--------------------------------------------------------------    
                [u_du,f_du,~,output] = fmincon(func,u_0,[],[],[],[],lb_local,ub_local, func_constr, options,problem_fix_d.par_objfun); %unconstrained
                
                if isempty(output.constrviolation)
                    c_du = 0;
                else
                    c_du = output.constrviolation;
                end
                
                if isempty(f_du)
                    f_du = f_d;
                end

                
                nfeval = nfeval + output.funcCount;
            else
                u_du = u_0;
                
                
                    try 
                      [f_du, c_du, ceq_du] = func(u_du, problem_fix_d.par_objfun);
                      objconstr = 1;
                    catch
                      f_du = func (u_du, problem_fix_d.par_objfun);
                    end
                
                %----------------------------------------------------------
                % CONSTRAINTS
                if ~isempty(problem_fix_d.constraint)                 % problem_fix_d.par_objfun.constraints_flag  == 1                   
                    
                    c_du = feval(problem_fix_d.constraint, u_du, problem_fix_d.par_objfun); % @mask_constraints_macsminmax_inner %[g,h]
                                 
                end
                %----------------------------------------------------------
                
                
                nfeval = nfeval + 1;
            end
            
            
            if f_du < f_d  &&  (c_du <= c_archive(idx_d-1) || c_du <=0)*(~isempty(problem_fix_d.constraint))
                f_d = f_du;
                u_d = u_du;
            end
            if (nargout > 3)
                all_f{obj}(idx_d,idx_u) = -minmaxsign*f_du;
            end
        end

        ARCHIVE_crosscheck{1, 1}{idx_d, 4} = -minmaxsign*f_d;
        ARCHIVE_crosscheck{1, 1}{idx_d, 2} =  map_affine(u_d, map_info);
%         f(idx_d,obj) = -minmaxsign*f_d;
%         u{obj}(idx_d,:) = u_d;
    end
end


return