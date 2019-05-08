function [ status ] = printSolutions( solutions, numb )
% [status ] = printSolutions( solutions, numb )
% 
% INPUTS:
%   - solutions Solution Dataset
%   - numb      Number of solutions to be displayed 
%
%--------------------------------------------------------------------------

if nargin > 1
    num_solutions = numb; 
else
    num_solutions = length(solutions);
end

% Print the Header
printHeader();

% Go throughout all the solutions
for i = 1 : num_solutions
   
    % Print the Current Solution
    disp([ sprintf('%6d | ', i) printSequence(solutions{i}, false)]);
   
    if i ~= num_solutions 
        disp('----------------------------------------------------------------------')
    end 
end

disp('======================================================================')