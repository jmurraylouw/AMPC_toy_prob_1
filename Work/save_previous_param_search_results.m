%% NB, must be for same dataset and same sigma
m=2;
% % 
% N_train_saved = 1100:100:1300; % Previously saved list of N to compute model for
% MAE_saved = rand(m,length(N_train_saved)); % Save best MAE for each N
% p_saved = randi(40,1,length(N_train_saved)); % Save p of best MAE for each N
% q_saved = randi(500,1,length(N_train_saved)); % Save q of best MAE for each N
% time_saved = rand(1,length(N_train_saved)); % Save time taken for each N
% 
N_train_list = 50:150:2600; % List of N to compute model for
MAE_list = rand(m,length(N_train_list)); % Save best MAE for each N
p_list = randi(40,1,length(N_train_list)); % Save p of best MAE for each N
q_list = randi(500,1,length(N_train_list)); % Save q of best MAE for each N
time_list = rand(1,length(N_train_list)); % Save time taken for each N



%% Merge new and saved results
% If two results are for the same N_train, keep the one with the smallest error

load('Data\N_train_error_time_HAVOK_sig=0.mat');

saved = [N_train_saved', mean(MAE_saved',2)]
list = [N_train_list', mean(MAE_list', 2)]

li = 1; % List index
si = 1; % Saved index

for k=1:1e10 % Basically while true
    if N_train_list(li) < N_train_saved(si) % list smaller
        N_train_saved   = insert(N_train_saved, si, N_train_list(li));
        MAE_saved       = insert(MAE_saved, si, MAE_list(:,li));
        p_saved         = insert(p_saved, si, p_list(li));
        q_saved         = insert(q_saved, si, q_list(li));
        time_saved      = insert(time_saved, si, time_list(li));
 
        si = si+1; % increment indices
        li = li+1;
    
    elseif N_train_list(li) == N_train_saved(si) % list equal
        if mean(MAE_list(:,li)) < mean(MAE_saved(:,si)) % If list has smaller MAE
            MAE_saved(:,si) = MAE_list(:,li);
            p_saved(si) = p_list(li);
            q_saved(si) = q_list(li);
            time_saved(si) = time_list(li);
        end
        
        si = si+1; % Increment indices
        li = li+1;
    
    else % list bigger
        si = si+1; % Keep comparing current li, only increment si
    end
    
    s_max = length(N_train_saved);
    l_max = length(N_train_list);
    
    if li > l_max % If all items in list checked
        break; 
    
    elseif si > s_max % If all in saved compared already
        N_train_saved = [N_train_saved, N_train_list(li:end)]; % Append rest of list to saved
        MAE_saved = [MAE_saved, MAE_list(:, li:end)]; % Append rest of list to saved
        p_saved = [p_saved, p_list(li:end)]; 
        q_saved = [q_saved, q_list(li:end)]; % Append rest of list to saved
        time_saved = [time_saved, time_list(li:end)]; % Append rest of list to saved   
        break;
    end
    
end

new_saved = [N_train_saved', mean(MAE_saved', 2)]

% Save the results to append to later
save('Data\N_train_error_time_HAVOK_sig=0.mat', 'N_train_saved', 'MAE_saved', 'p_saved', 'q_saved', 'time_saved');

bar(mean(MAE_saved',2))

function new_array = insert(array, index, entry)
    if index == 1 % to avoid index-1 = 0
        new_array = [entry, array];
    else
        new_array = [array(:, 1:index-1), entry, array(:, index:end)];
    end
end

