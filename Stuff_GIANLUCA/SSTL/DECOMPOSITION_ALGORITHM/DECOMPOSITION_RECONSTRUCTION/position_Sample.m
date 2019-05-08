function [PS] = position_Sample(num_sample_tot, num_coupled_vectors, num_sample, Sample, in)


% [i(1),i(2),i(3),i(4),i(5),i(6),i(7),i(8)] = ind2sub(size(max_mat), num_sample);

A =1;
for num_Belief = num_coupled_vectors:-1:2
    if in.dim_u(in.num_functions + num_Belief)>0
    if in.output == 0 || in.output == 2
        A = A*length(Sample{num_Belief}.FE);
    elseif in.output == 1
        A = A*length(Sample{num_Belief}.FE_Plausibility);
    end
    PS(num_Belief) = ceil(num_sample*A/num_sample_tot);
    num_sample = num_sample-num_sample_tot/A*(PS(num_Belief)-1);
    end
end
PS(1) = num_sample;



end