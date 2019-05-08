function [ solutions_out best_seq] = sortSolution( solutions_in, wa, wb)
%  SortSolution(solutions, wa, wb)
%
% INPUTS:
%   
%  - solutions The solution to be sorted
%  - wa  Weight for the total deltaV
%  - wb  Weight for the total lenght of the sequence  
%
% OUTPUTs:
%
%  solutions_out The sorted solutions
%  best_seq      The best solution acording the given weights
%
%--------------------------------------------------------------------------


solutions_out = {};

% Sanity Check
if length(solutions_in) <= 0
    warning('The Solutions_in is empty.')
    return; 
end

% Sanity Check
if ~isa(solutions_in{1}{end}, 'TrajectoryArc')      
    warning('The Solutions_in is not storing TrajectoryArc Elements.')
    return; 
end


% Find the Maximum deltaV
maxdV       = -inf;
maxElements = -inf;
best_seq = {};
temp_totaldV = inf;

% Search for the maximum total delta-V and the longer sequence
for i = 1:  length(solutions_in)
    
    % Get Current total deltaV
     curr_dV = solutions_in{i}{end}.total_delatV;

     if curr_dV > maxdV
         maxdV = curr_dV;
     end
     
     % Get Current total deltaV
     curr_numbElements = solutions_in{i}{end}.id;

     if curr_numbElements > maxElements
         maxElements = curr_numbElements;
     end
          
     if length(solutions_in{i}) == 6           
         if ( solutions_in{i}{end}.total_delatV < temp_totaldV )
             best_seq = solutions_in{i};
             temp_totaldV = solutions_in{i}{end}.total_delatV;
         end
     end     
end


% Pre-allocate the Weights Array
weights = ones(1,length(solutions_in));

% Populate the Weights Array
for i = 1:  length(solutions_in)
    
    % Get the Current Total deltaV
    curr_totaldV = solutions_in{i}{end}.total_delatV/maxdV;
    
    % Get the Current number of visited Asteroids
    curr_numbEle = solutions_in{i}{end}.id/maxElements;
   
    % Compute the Curent Weight Value
    weights(i) =  wa*curr_totaldV + wb*(1/curr_numbEle);
    
end

% Sort the Weight Array
[temp, indixes] = sort(weights);

% Translate the sorted weight array to the solutions
%solutions_out = solutions_in(fliplr(indixes));
solutions_out = solutions_in((indixes));

