function [sample] = presample_fe(n,n_int)
% returns a sample of cardinality n in the focal element structure defined by n_int
% n_int is a vector [n1, n2, ...] where ni is the number of divisions in dimension i
% coded very fast so it's not good when n approaches cardinality of domain


ndim = length(n_int);
buckets = cell(1,ndim);
sample = [];
m = n;
i0 = 1;
while m>0
    for i = i0:n
        for j = 1:ndim
            if (isempty(buckets{j}))
                buckets{j} = randperm(n_int(j));
            end
            sample(i,j) = buckets{j}(end);
            buckets{j} = buckets{j}(1:end-1);
        end
    end
    sample = unique(sample,'rows');
    i0 = size(sample,1)+1;
    m = n - i0 + 1;
end

return
