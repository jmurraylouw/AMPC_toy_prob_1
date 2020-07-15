disp('-------------------------')
%% Load records
try
    load(search_space_file);
catch
    disp('No search space file')
    N_train_record = [];
    q_record = [];
    p_record = [];    
    search_space = zeros(0, 0, 0);
end

% N_train_record = 100:100:200;
% q_record = 10:12;
% p_record = 8:10;
% 
% search_space = zeros(length(N_train_record), length(q_record), length(p_record))

%% Search lists
N_train_list = 200:100:400;
q_search = [1];
p_search = [6,5,8];

%% for N_train
for N_train = N_train_list
    N_train_i = find(N_train_record == N_train); % index of N_train in search_space
    if isempty(N_train_i) % If N_train not in record list
        N_train_i = length(N_train_record)+1; % q index becomes next in list
    end
    
    %% for q
    for q = q_search
        q_i = find(q_record == q); % index of q in search_space
        if isempty(q_i) % If q not in record list
            q_i = length(q_record)+1; % q index becomes next in list
        end
        
        %% for p
        for p = p_search
            p_i = find(p_record == p); % index of p in search_space
            if isempty(p_i) % If p not in record list
                p_i = length(p_record)+1; % q index becomes next in list
            end
            
            %% Determine if new_search or not
            try
                MAE_record = search_space(N_train_i, q_i, p_i);
                if MAE_record == 0 % Zero error means unsearched
                    new_search = 1; % 1 = first time searching this param combo
                else
                    new_search = 0; % 0 = has been searched before
                end
            catch % Catch index out of bounds
                new_search = 1; % 1 = first time searching this param combo
            end
            
            if new_search
                %% Do p calculations...
                
                %% Record search
                search_space(N_train_i, q_i, p_i) = MAE;
                N_train_record(N_train_i) = N_train; % record q in search
                q_record(q_i) = q; % record q in search
                p_record(p_i) = p; % record q in search
            else
                % Do nothing   
            end            
        end % p
    end % q
end % N_train

%% Save search_space
search_space
N_train_record
q_record
p_record
save(search_space_file, 'search_space', 'N_train_record', 'q_record', 'p_record')






