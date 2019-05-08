%==========================================================================
% This Script Generate a graph showing the evolution of the DeltaV of the
% solutions already organized in groups. 
%==========================================================================

% Set the Solutions
groups = {solutions_sorted4, solutions_sorted5, solutions_sorted6};

% 
deltaVs = {};

% (1) Go through out all the groups
for g = 1: length(groups)
    
    % Get the current Group
    curr_group = groups{g};
    
    curr_deltaVs = [];
    
    % (2) Go through out all the sequences on the group
    for i = 1: length(curr_group)
        
        % Get the Current Sequence
        curr_seq = curr_group{i};
        
        % Store the deltaV for the current Sequence
        curr_deltaVs(i) = curr_seq{end}.total_delatV;
                         
    end % End (2)
    
    deltaVs{end + 1} = curr_deltaVs;
    
end % End (1)
% 
% % Plot the solutions
% plot(deltaVs{1, 1})
% hold on
% plot(deltaVs{1, 2})
% 
% 
% x = [1:length(deltaVs{1, 1})];
% hold on
% plot(x, deltaVs{1, 1}(1));
% hold on
% plot(x, deltaVs{1, 1}(end));
% hold on
% plot(x, deltaVs{1, 2}(1));
% hold on
% plot(x, deltaVs{1, 2}(end));


hold off

[n,xout] = hist(deltaVs{1,1},20);
bar(xout,n)

[n,xout] = hist(deltaVs{1,2},20);
hold on
bar(xout,n)

[n,xout] = hist(deltaVs{1,3},20);
hold on
bar(xout,n)

