function index=cartprod(n)
% Cartesian product matrix.
%
%   index = cartprod(n)
%
% INPUT
%   n       Integer vector with the number of indexes for each variable.
%
% OUTPUT
%   index   Cartesian product matrix. One element for each row. It has a
%           number of columns equal to length(n), and a number of rows
%           equal to prod(n).
%
% EXAMPLE
%   index = cartprod([2 3 2]) gives:
%       index = [1     1     1
%                1     1     2
%                1     2     1
%                1     2     2
%                1     3     1
%                1     3     2
%                2     1     1
%                2     1     2
%                2     2     1
%                2     2     2
%                2     3     1
%                2     3     2]
%
% FUNCTIONS CALLED
%   (none)
%
% Matteo Ceriotti, 15-06-2006
% - Revised by Camilla Colombo - 14-05-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

npar=length(n);
nindex=1;
for i=1:npar
    nindex=nindex*n(i);
end
index=ones(nindex,npar);
index(1,npar)=0;
indexcount=0;
for i=1:nindex
    index(i,npar)=indexcount+1;
    j=npar;
    while (j>0)
        if(index(i,j)>n(j))
            index(i:nindex,j)=1;
            index(i:nindex,j-1)=index(i:nindex,j-1)+1;
        end
        j=j-1;
    end
    indexcount=index(i,npar);
end

return