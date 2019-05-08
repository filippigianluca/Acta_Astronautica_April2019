function [out] = save_MP_AIDEA_closefile(algo_inner)


    % SAVE POPULATION MP_AIDEA  inner
    
    if algo_inner.par.save_pop_DE
        for iii = 1 : algo_inner.par.n_populations
            fclose(algo_inner.par.fileID(iii))
        end
    end
    
    if algo_inner.par.save_local_search
        for iii = 1 : algo_inner.par.n_populations
            fclose(algo_inner.par.fileID2(iii))
        end
    end
    
    if algo_inner.par.save_pop_LR
        for iii = 1 : algo_inner.par.n_populations
            fclose(algo_inner.par.fileID3(iii))
        end
    end
    if algo_inner.par.save_pop_GR
        for iii = 1 : algo_inner.par.n_populations
            fclose(algo_inner.par.fileID4(iii))
        end
    end

out = 1;
return