% Clean the Command Window
clc

% Sorting Flag
flag_sort = false;

% Weight Factor for the total deltaV
wdv = 1; 

% Weight Factor for the Sequence Length
wse = 0;

% Number solution to Print
num_solution = 5000;

% Sorte the Solutions 
if flag_sort
    solutions_sorted6 = sortSolution(groups_6, wdv, wse);
end

% Display the solutions
printSolutions(solutions_sorted6, num_solution);