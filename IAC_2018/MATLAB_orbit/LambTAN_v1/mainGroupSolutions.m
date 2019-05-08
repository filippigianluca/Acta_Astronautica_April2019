% Create Histogram Analysis. 
%
% This function run through a set of sequences (solutions) and create a
% histogram of the number of visit asteroids

clear groups_4 groups_5 groups_6

groups_4 = {};
groups_5 = {};
groups_6 = {};

% Go through out all the solutions
for i = 1 : length(solutions)   
    
    % Get the Current Sequence
    sequence = solutions{i};
    
    % Create the key 
    key = length(sequence);
    
    switch key 
        case 6 
            groups_6{end + 1} = sequence;
        case 5 
            groups_5{end + 1} = sequence;
        case 4
            groups_4{end + 1} = sequence;
    end         
end

%
solutions_sorted6 = sortSolution(groups_6, 1.0, 0.0);
solutions_sorted5 = sortSolution(groups_5, 1.0, 0.0);
solutions_sorted4 = sortSolution(groups_4, 1.0, 0.0);

%--------------------------------------------------------------------------
% Step 2: 
%--------------------------------------------------------------------------

if(length(groups_6) > 0) 
    
    listSol = containers.Map;
    
    short_solutions_sorted6 = {};
    
    % Go through out all the solutions
    for i = 1 : length(solutions_sorted6)
        
        % Get the Current Sequence
        sequence = solutions_sorted6{i};
        
        % Create the key
        key_sequence = 'Earth ';
        
        for i = 1 : length(sequence)
            key_sequence = [key_sequence '  ' sequence{i}.target_celbod.name];
        end
        
        % check if the sequence exists
        if listSol.isKey(key_sequence)
            continue;
        end
        
        listSol(key_sequence) = sequence;
        short_solutions_sorted6{end + 1} = sequence;
    end
    
    printSolutions(solutions_sorted6, 20);
    printSolutions(short_solutions_sorted6);
end 
%--------------------------------------------------------------------------
% Step : 
%--------------------------------------------------------------------------
listSol2 = containers.Map;

short_solutions_sorted5 = {};

% Go through out all the solutions
for i = 1 : length(solutions_sorted5)   
    
    % Get the Current Sequence
    sequence = solutions_sorted5{i};
    
    % Create the key 
    key_sequence = 'Earth ';
 
    for i = 1 : length(sequence)
        key_sequence = [key_sequence '  ' sequence{i}.target_celbod.name];
    end
    
    % check if the sequence exists
    if listSol2.isKey(key_sequence)
        continue;
    end
    
    listSol2(key_sequence) = sequence;
    short_solutions_sorted5{end + 1} = sequence;
end

printSolutions(solutions_sorted5, 20);
printSolutions(short_solutions_sorted5);

%--------------------------------------------------------------------------
% Step : 
%--------------------------------------------------------------------------
listSol3 = containers.Map;

short_solutions_sorted4 = {};

% Go through out all the solutions
for i = 1 : length(solutions_sorted4)   
    
    % Get the Current Sequence
    sequence = solutions_sorted4{i};
    
    % Create the key 
    key_sequence = 'Earth ';
 
    for i = 1 : length(sequence)
        key_sequence = [key_sequence '  ' sequence{i}.target_celbod.name];
    end
    
    % check if the sequence exists
    if listSol3.isKey(key_sequence)
        continue;
    end
    
    listSol3(key_sequence) = sequence;
    short_solutions_sorted4{end + 1} = sequence;
end


printSolutions(solutions_sorted4, 20);
printSolutions(short_solutions_sorted4);
