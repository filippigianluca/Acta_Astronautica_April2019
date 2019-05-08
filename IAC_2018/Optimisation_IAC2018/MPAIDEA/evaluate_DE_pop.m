function [Val, ibest, BestVal, nFeVal, exit_DE] = ...
    evaluate_DE_pop(i_pop_number, pop, fname, nFeVal, nFeValMax, ibest, BestVal, Val)



exit_DE.flag = 0;

% Call runobjconstr for global part of the algorithm
flag_LG.global = 1;
flag_LG.local  = 0;

% The population is composed by NP 
[NP,~]    = size(pop);



nFeVal_temp = nFeVal(1,i_pop_number);



%%%%% PARFOR
% Problem: I cannot check if I have reached the maximum number of function
% evaluations during the parfor
parfor i = 2 : NP                        % check the remaining members
    
    
    % Evaluation of the function: different outputs depending on the presence
    % of constraints and on how the constraints are handled (see above for
    % additional comments)
    if ( fname.weighted && ~isempty(fname.constr) ) || isempty(fname.constr)
        [~,Val(i)] = runobjconstr(pop(i,:), fname, flag_LG, [], [], []);
    elseif fname.weighted == 0 && ~isempty(fname.constr)
        [~,Val_temp(i)] = runobjconstr(pop(i,:), fname, flag_LG, [], [], []);
    end
    
    nFeVal_temp = nFeVal_temp + 1;
    

end
% =====================================================================

% =====================================================================
% If maximum number of function evaluation is reached
if sum(nFeVal_temp) >= nFeValMax
    exit_DE.new_elements = NP;
    exit_DE.iter = 0;
    
    % -----------------------------------------------------------------
    % If constraints are not weighted, then it is necessary to find the maximum
    % value of the objective function
    if fname.weighted == 0 && ~isempty(fname.constr)
        
        % If there are feasible individuals...
        if any( ( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
                ( cell2mat({Val_temp.non_feas_c})  == 0 ) )
            % Feasible individuals: objective function is objective function
            % Feasible individuals are denoted by
            % ( cell2mat({Val_temp.non_feas})  == 0 )
            % thet is, the flag for infeasibility if put to zero
            Val( logical (( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
                ( cell2mat({Val_temp.non_feas_c})  == 0 ) ) ) = ...
                cellfun(@max, {Val_temp(logical (( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
                ( cell2mat({Val_temp.non_feas_c})  == 0 ) )).yy})';
            % Above, cellfun(@max, {}) is used because otherwise only
            % the last value is assigned to all the Val - bug found by
            % Gianluca
            
        end
        
        % If there are unfeasible individuals, either for ceq or c..
        if any( ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
                ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0 )
            % Non feasible individuals: the objective function is the maximum of
            % the objective function for all the individuals + equality
            % constraint + maximum of inequality constraints for that individual
            % Non feasible individuals are denoted by index
            % ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
            %    ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0
            % thet is, the flag for infeasibility if put to 1
            index_nonfeas = logical( ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
                ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0 );
            % Objective: maximum of objective function + violation of
            % equality constraint (only if violated) + violation of
            % inequality constraint (only if violated)
            Val( index_nonfeas ) = ...
                max(cell2mat({Val_temp.yy})) + ...
                cellfun(@norm,{Val_temp( index_nonfeas ).ceq})' .* ...
                cellfun(@norm, {Val_temp( index_nonfeas ).non_feas_ceq})' + ...
                cellfun(@max, {Val_temp( index_nonfeas ).c})' .* ...
                cellfun(@norm, {Val_temp( index_nonfeas ).non_feas_c})';
        end
        
        [BestVal, ibest] = min(Val);
        
    end
    % -----------------------------------------------------------------
    
    %         BestMem = pop(ibest,:);
    
    exit_DE = 1;
end
% end


% if member is better
if ( fname.weighted || isempty(fname.constr) )
    [BestVal_temp, ibest_temp] = min(Val);
    if BestVal_temp < BestVal
        BestVal = BestVal_temp;
        ibest = ibest_temp;
    end
    %         BestMem = pop(i,:);
end

% Now all the objectives values are available for all the individuals...
% If constraints are not weighted, then it is necessary to find the maximum
% value of the objective function
if fname.weighted == 0 && ~isempty(fname.constr)  
    
    % If there are feasible individuals (both equality and inequality feasibility)...
    if any( ( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
            ( cell2mat({Val_temp.non_feas_c})  == 0 )) 
        % Feasible individuals: objective function is objective function
        % Feasible individuals are denoted by
        % ( cell2mat({Val_temp.non_feas})  == 0 ) .* ...
        %     ( cell2mat({Val_temp.non_feas_c})  == 0 ) == 1
        % thet is, the flag for infeasibility if put to zero for both
        % equality and inequality
        Val(  logical (( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
            ( cell2mat({Val_temp.non_feas_c})  == 0 ) ) ) = ...
            cellfun(@max, {Val_temp( logical (( cell2mat({Val_temp.non_feas_ceq})  == 0 ) .* ...
            ( cell2mat({Val_temp.non_feas_c})  == 0 ) ) ).yy})';
        % In order to take all the variables, cellfun(@max, {}) is used -
        % bug found by Gianluca
        
    end
    
    % If there are unfeasible individuals..
    if any( ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
                    ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0 )
        % Non feasible individuals: the objective function is the maximum of
        % the objective function for all the individuals + equality
        % constraint + maximum of inequality constraints for that individual
        % Non feasible individuals are denoted by index
        % ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
                %    ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0 
        % thet is, the flag for infeasibility if put to zero
        index_nonfeas = logical( ( cell2mat({Val_temp.non_feas_ceq})  ~= 0 ) + ...
                    ( cell2mat({Val_temp.non_feas_c})  ~= 0 ) ~=0 );
        % Objective: maximum of objective function + violation of
        % equality constraint (only if violated) + violation of
        % inequality constraint (only if violated)
        Val( index_nonfeas ) = ...
            max(cell2mat({Val_temp.yy})) + ...
            cellfun(@norm,{Val_temp( index_nonfeas ).ceq})' .* ...
            cellfun(@norm, {Val_temp( index_nonfeas ).non_feas_ceq})' + ...
            cellfun(@max, {Val_temp( index_nonfeas ).c})' .* ...
            cellfun(@norm, {Val_temp( index_nonfeas ).non_feas_c})';
    end
    
    [BestVal, ibest] = min(Val);
%     BestMem = pop(ibest,:);
    
end


nFeVal(1,i_pop_number) = nFeVal_temp;