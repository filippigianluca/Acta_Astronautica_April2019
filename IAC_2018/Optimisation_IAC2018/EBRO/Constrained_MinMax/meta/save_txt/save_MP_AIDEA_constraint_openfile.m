function algo_inner = save_MP_AIDEA_constraint_openfile(problem_minmax, N_u_mpaidea, algo_inner)



% SAVE POPULATION MP_AIDEA  constraints

algo_inner.par.save_pop_DE = algo_inner.con.save_pop_DE;
% If yes, choose a name for the file:
name_save_pop_DE =['population_DE_algo_constraint' num2str(N_u_mpaidea) '.txt'];


algo_inner.par.save_local_search = algo_inner.con.save_local_search;
% If yes, choose name for file:
name_save_local_search = ['minima_fmincon_constraint' num2str(N_u_mpaidea) '.txt'];


% -------------------------------------------------------------------------
% Save populations at local restart (each one saved on a different file)?
% -------------------------------------------------------------------------
algo_inner.par.save_pop_LR = algo_inner.con.save_pop_LR;
% If yes, choose prefix of name for files:
name_save_pop_LR = ['pop_LR_constraint' num2str(N_u_mpaidea)];

% -------------------------------------------------------------------------
% Save populations at global restart (each one saved on a different file)?
% -------------------------------------------------------------------------
algo_inner.par.save_pop_GR = algo_inner.con.save_pop_GR;
% If yes, choose prefix of name for files:
name_save_pop_GR = ['pop_GR_constraint' num2str(N_u_mpaidea)];



% File to save population of DE
if algo_inner.par.save_pop_DE
    algo_inner.par.str = '%8.6e';
    for iii = 1 : problem_minmax.dim_u
        algo_inner.par.str = [algo_inner.par.str,' ', '%8.6e'];
    end
    algo_inner.par.str = [algo_inner.par.str, '\n'];
        for s = 1 : algo_inner.par.n_populations
            algo_inner.par.fileID(s) = fopen(strcat(name_save_pop_DE, num2str(s),'.txt'),'w');
        end     
%         algo_inner.par.fileID = fopen(name_save_pop_DE,'w');
end


% File to save local searches
if algo_inner.par.save_local_search
    % If yes, uncomment the following and give name to files:
    algo_inner.par.str = '%8.6e';
    for iii = 1 : problem_minmax.dim_u
        algo_inner.par.str = [algo_inner.par.str,' ', '%8.6e'];
    end
    algo_inner.par.str = [algo_inner.par.str, '\n'];
    for s = 1 : algo_inner.par.n_populations
        algo_inner.par.fileID2(s) = fopen(strcat(name_save_local_search, num2str(s),'.txt'),'w');
    end    
%         algo_inner.par.fileID2 = fopen(name_save_local_search,'w');
end

% File to save local restarts
if algo_inner.par.save_pop_LR
    % If yes, uncomment the following and give name to files:
    algo_inner.par.str2 = '%8.6e';
    for iii = 1 : problem_minmax.dim_u - 1
        algo_inner.par.str2 = [algo_inner.par.str2,' ', '%8.6e'];
    end
    algo_inner.par.str2 = [algo_inner.par.str2, '\n'];
    for iii = 1 : algo_inner.par.n_populations
        algo_inner.par.fileID3(iii) = fopen(strcat(name_save_pop_LR,num2str(iii),'.txt'),'w');
    end
end

% File to save global restarts
if algo_inner.par.save_pop_GR
    % If yes, uncomment the following and give name to files:
    algo_inner.par.str2 = '%8.6e';
    for iii = 1 : problem_minmax.dim_u - 1
        algo_inner.par.str2 = [algo_inner.par.str2,' ', '%8.6e'];
    end
    algo_inner.par.str2 = [algo_inner.par.str2, '\n'];
    for iii = 1 : algo_inner.par.n_populations
        algo_inner.par.fileID4(iii) = fopen(strcat(name_save_pop_GR,num2str(iii),'.txt'),'w');
    end
end


return
