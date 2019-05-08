function [str_archive, str_archive_inner, str_archive_outer] = open_txt_optimise_so(problem_minmax, algo_inner)


    str_archive       = [];
    str_archive_inner = [];
    str_archive_outer = [];
    
    
    %----------------------------------------------------------------------
    % TXT FILES
    %----------------------------------------------------------------------    
    if algo_inner.save_archive
        %------------------------------------------------------------------
        fileID_archive = fopen('ARCHIVE.txt','w');
        str_archive = [];
        for len_arch = 1: sum([problem_minmax.dim_u   problem_minmax.dim_d   3])
            str_archive = [str_archive,' ', '%18.16f'];
        end
        str_archive = [str_archive, '\n'];
        
        %------------------------------------------------------------------
        fileID_archive_inner = fopen('ARCHIVE_INNER.txt','w');
        str_archive_inner = [];
        for len_arch = 1:2
            str_archive_inner = [str_archive_inner,' ', '%18.16f'];
        end
        str_archive_inner = [str_archive_inner, '\n'];
        
        %------------------------------------------------------------------
        fileID_archive_outer = fopen('ARCHIVE_OUTER.txt','w');
        str_archive_outer = [];
        for len_arch = 1:2
            str_archive_outer = [str_archive_outer,' ', '%18.16f'];
        end
        str_archive_outer = [str_archive_outer, '\n'];
    end
    %----------------------------------------------------------------------
    
end