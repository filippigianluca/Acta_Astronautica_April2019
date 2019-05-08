function out = cross_fast(r1,r2)
% Very fast cross product. It requires about 1/10 of the time of cross.m.
% The output is a row or column vector, as in the first input vector.
%
%   out = cross_fast(r1, r2)
%
% INPUT
%   r1      First vector.
%   r2      Second vector.
%
% OUTPUT
%   out     Cross product between r1 and r2. Row or column depending on r1.
%
% FUNCTIONS CALLED
%   (none)
%
% Matteo Ceriotti, 12-02-2007
% Revised by Camilla Colombo, 14-02-2007

if size(r1,1)==1
    out=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)];
else
    out=[r1(2)*r2(3)-r1(3)*r2(2),r1(3)*r2(1)-r1(1)*r2(3),r1(1)*r2(2)-r1(2)*r2(1)]';
end