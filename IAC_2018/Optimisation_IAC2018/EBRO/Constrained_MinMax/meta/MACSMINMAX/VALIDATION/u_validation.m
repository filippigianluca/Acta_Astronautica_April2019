% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [output_validation,u,nfeval,all_f] = u_validation(problem_fix_d,d_record,u_record,local_search_flag,objectives)

objconstr = 0;

minmaxsign = problem_fix_d.par_objfun.sign;
archsize = size(d_record, 1);
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

for idx_d = 1:archsize
    % fix d and scale
    problem_fix_d.par_objfun.d = d_record(idx_d,:);
    for obj = objectives
        problem_fix_d.par_objfun.objective = obj;
        f_du = [];

        u_du=[];
        c_d  = nan;
        f_d=nan;
        u_d=nan(1,problem_fix_d.dim);
        for idx_u = 1:size(u_record{obj},1)
            u_0=u_record{obj}(idx_u,:);
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
                if ~isempty(problem_fix_d.constraint)             
                    % CONSTRAINT
%                     [constraint] = feval(problem_fix_d.par_objfun.mask_constraints, u_du, problem_fix_d.par_objfun);
                    func_constr = problem_fix_d.par_objfun.mask_constraint;
                else
                    func_constr = [];
                end
                %--------------------------------------------------------------    
                [u_du,f_du,~,output] = fmincon(func,u_0,[],[],[],[],lb_local,ub_local,func_constr,options,problem_fix_d.par_objfun); %unconstrained
                    
                constraint = output.constrviolation;
                
                nfeval = nfeval + output.funcCount;

                    
                
                
            else
                
%                 if problem_fix_d.fitnessfcn.obj_constr == 0 % function and constraint are in different functions
                    u_du = u_0;
                    
                    try 
                      [f_du, c_du, ceq_du] = func(u_du, problem_fix_d.par_objfun);
                      objconstr = 1;
                    catch
                      f_du = func(u_du, problem_fix_d.par_objfun);
                    end
                    
                    
%                     f_Au = func(u_du, problem_fix_d.par_objfun);                    
%                     if isstruct(f_Au)
%                         f_du = f_Au.f;
%                         c_du = f_Au.c;
%                         ceq_du = f_Au.ce;
%                     else
%                         f_du = f_Au;
%                     end
                    
                    
                    
%                     f_du = func(u_du, problem_fix_d.par_objfun);
                    %----------------------------------------------------------
                    % CONSTRAINTS
%                     if ~isempty(problem_fix_d.constraint)
%                         [constraint] = feval(problem_fix_d.constraint, u_du, problem_fix_d.par_objfun);
%                     end
                    %----------------------------------------------------------
                    nfeval = nfeval + 1;
                
%                 else % function and constraint are in the same function
%                     u_du = u_0;
%                     [f_du, constraint, ~] = func(u_du, problem_fix_d.par_objfun);
%                     nfeval = nfeval + 1;                   
%                 end
            end
            
            
            
            if objconstr
                if (c_du > c_d) || (idx_u==1)    % -constraint>=0&&
                    c_d = c_du;
                end          
            end
            
            
            
            
%             if isempty(problem_fix_d.constraint)
                if (f_du < f_d) || (idx_u==1)    % -constraint>=0&&
                    f_d = f_du;
                    u_d = u_du;
                end
            
            %--------------------------------------------------------------
            % CONSTRAINT
            %
            % if the constraint is violated
%             elseif ~isempty(problem_fix_d.constraint)
%                 
%                 constraint_i(idx_u) = constraint;
%                 f_du_i(idx_u) = f_du;
%                 if constraint > 0
                %                 if constraint > 0 && (f_du - constraint < f_d)
%                 if constraint > 0 && - constraint < f_du

% % % % % %                 if any(constraint_i <= 0)
% % % % % %                     [~, b] = find(constraint_i <= 0);
% % % % % % 
% % % % % %                     [f_d, c] = min(f_du_i(b));
% % % % % %                     u_d = u_record{obj}(c,:);
% % % % % %                 else
% % % % % %                     [~, b] = min(constraint_i);
% % % % % %                     f_d = f_du_i(b);
% % % % % %                     u_d = u_record{obj}(b,:);
% % % % % %                 end
    
%                 if (f_du < f_d) || (idx_u==1)    % -constraint>=0&&
%                     f_d = f_du;
%                     u_d = u_du;
%                 end
%                     f_d = -max(0, constraint) + min(f_du, f_d);
%                     u_d = u_du;
%                 end
%                 end
%             end
            %--------------------------------------------------------------
            
            if (nargout > 3)
                all_f{obj}(idx_d,idx_u) = -minmaxsign*f_du;
            end
            
            
            
        end
        
        
        u{obj}(idx_d,:) = u_d;

if objconstr 
    output_validation.f(idx_d,obj) = -minmaxsign*f_d;
    output_validation.c(idx_d,obj) = c_d;
    output_validation.ceq = [];
else
    output_validation(idx_d,obj) = -minmaxsign*f_d;
end
        

    end
end

return