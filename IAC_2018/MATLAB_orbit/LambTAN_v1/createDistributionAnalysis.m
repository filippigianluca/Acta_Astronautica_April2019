function [] = createDistributionAnalysis(solutions)
% Create Distribution Analysis. 
%
% This function compute the Distribution of the visited asteroids within all 
% the solutions.
%
% PROTOTYPE:
%
%     createDistributionAnalysis(solutions)
%
% INPUTS:
%    solutions  The Solution Cell Array
%
% OUTPUT:
%
% No values are returned. The results are displayed in the Command Window
% with the following format:
%
%
% Asteroid NAME_OF_ASTEROID _________________________________
%      A:   B% (NUM)	
%
% Where:
%      A    The position on the sequence
%      B    Percentage of how much this current asteroid has found in all 
%           the sequences 
%      NUM  Total number of astroid found in this sequence position
%
% -------------------------------------------------------------------------

% Create the HashMap
hmapAsteroids = containers.Map();

% Go through out all the solutions
for i = 1 : length(solutions)   
    
    % Get the Current Sequence
    sequence = solutions{i};
    
    % Go through out the sequence
    for j = 1 : length(sequence)
                
        % Get the Name of Asteroid
        asteroid = sequence{j}.target_celbod.name;
                
        % Increase the counter in the HashMap
        if hmapAsteroids.isKey(asteroid)
            
            hmap = hmapAsteroids(asteroid);
            key  = sprintf('%d',j);
            
            if hmap.isKey(key)
                hmap(key) = hmap(key) + 1;
            else
                hmap(key) = 1;       
            end
            
            hmap('total') = hmap('total') + 1;  
            
        else
            hmapAsteroids(asteroid) = containers.Map();
            hmap = hmapAsteroids(asteroid);
            hmap('total') = 1;            
        end
    end        
end
 
% % Plot Histogram
% bar( cell2mat(histogram.values));
% set(gca, 'XTick', 1:histogram.Count, 'XTickLabel', histogram.keys);

keys_asteroids = hmapAsteroids.keys;

% Go through out all the solutions
for i = 1 : length(keys_asteroids) 
   
    hmap2 = hmapAsteroids(keys_asteroids{i});
    disp(['Asteroid ' keys_asteroids{i} ' _________________________________']);
    
    curr_keys_asteroid = hmap2.keys;
    total = hmap2('total');
    sline = '';
    for j=1 : length(curr_keys_asteroid)
        curr_pos = hmap2(curr_keys_asteroid{j});
        sline = [sline sprintf(' %s: %d (%.1f%%)\t', curr_keys_asteroid{j}, curr_pos, curr_pos/total * 100 )];
        disp( sprintf(' %5s: %5.1f%% (%d)\t', curr_keys_asteroid{j}, curr_pos/total * 100, curr_pos ) );
    end
    
    %disp(sline);
    disp(' ');
end

end
