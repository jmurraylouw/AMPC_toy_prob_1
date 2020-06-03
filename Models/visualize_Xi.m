function vis_Xi = visualize_Xi(x_names, Xi, polyorder)
% Create cell array to visualise terms in Xi matrix from SINDY
% x_names   = string names of terms in Theta e.g {'x', 'y', 'z', 'sin(y)', '(1/x)'}
% Xi        = Xi matrix;
% polyorder = highest order of polynomial in function

num_terms = length(x_names);
vis_Xi = cell(size(Xi)+[1 2]); % Empty visualisation cell matrix
% Leave space in vis_Xi for 2 columns on left and one on top for headings

vis_Xi(2:end, 3:end) = num2cell(Xi); % Insert Xi into vis_Xi 
rows = size(Xi,1); % Number of rows in Xi/ states
cols = size(Xi,2); % Number of columns in Xi/ terms

for i = 1:cols % Add column headers of states
    vis_Xi(1,i+2) = strcat(x_names(i),'_dot');
end

vis_Xi(2,2) = {'1'};
index = 3; % Row index to add next label

for i=1:num_terms
    vis_Xi(index,2) = x_names(i);
    index = index+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:num_terms
        for j=i:num_terms
            vis_Xi{index,2} = [x_names{i},x_names{j}];
            index = index+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:num_terms
        for j=i:num_terms
            for k=j:num_terms
                vis_Xi{index,2} = [x_names{i},x_names{j},x_names{k}];
                index = index+1;
            end
        end
    end
end

% assert(index-2 == rows, 'There are %d rows in Xi, but only %d combinations of terms',rows,index-2 )
if(index-2 ~= rows) % Add denominator terms
    vis_Xi(((end-1)/2+1+1):end,2) = vis_Xi(2:((end-1)/2+1),2);
end

vis_Xi(:,1) = num2cell((0:size(vis_Xi,1)-1)'); % Add row numbering
