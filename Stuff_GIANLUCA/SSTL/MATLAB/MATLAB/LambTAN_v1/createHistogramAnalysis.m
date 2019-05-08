function [] = createHistogramAnalysis(solutions)
% Create Histogram Analysis. 
%
% This function run through a set of sequences (solutions) and create a
% histogram of the number of visit asteroids
%

% Create the HashMap
histogram = containers.Map();

% Go through out all the solutions
for i = 1 : length(solutions)   
    
    % Get the Current Sequence
    sequence = solutions{i};
    
    % Go through out the sequence
    for j = 1 : length(sequence)
                
        % Get the Name of Asteroid
        asteroid = sequence{j}.target_celbod.name;
                
        % Increase the counter in the HashMap
        if histogram.isKey(asteroid)
            histogram(asteroid) = histogram(asteroid) + 1;
        else
            histogram(asteroid) = 1;
        end
    end        
end

% Plot Histogram
bar( cell2mat(histogram.values));
set(gca, 'XTick', 1:histogram.Count, 'XTickLabel', histogram.keys);

end
