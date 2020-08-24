%% HAVOK N_train vs error sweep with grid search for parameters
% HAVOK with control - of cart pendulum
% Sweep through different N_train, 
% do grid search to find best hyper-parameters for each N_train with sigma
% save best error with corresponding hyper-parameters for each N_train
% save hyper-parameters already searched for each N_train and sigma
%   combination to avoid searching same search space again
% 
% Use corresponding HAVOK_script to load and run best hyperparameters to 
% see prediction results 

% clear all;
close all;
total_timer = tic;

%% Read data

% load('floating_pend_2D_random_1.mat') % Load simulation data
% load('floating_pend_2D_PI_z_control_and_random_1.mat') % simulation data with z controlled by PI controller
simulation_data_file = 'floating_pend_2D_data_1';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Extract data
u_data  = out.u.Data';
x_data  = out.x.Data';
y_data  = x_data([1:3],:); % Measurement data (x, z, theta)
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
u0 = [0; 6*9.81]; % Input needed to keep at a fixed point
u_data  = u_data - u0; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 5000; % Num of data samples for testing
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

%% Parameters
% Very dependant on choice of p, r, q

sigma = 0.01; % Noise standard deviation

state_scale = [1; 1; 1]; % Scale MAE of each state by this when comparing models

c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)

% N_train; % Num of data samples for training, rest for testing
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
% p = Truncated rank of system of Omega
% r = Truncated rank of system of X2
% q = number of delays

% Search parameters
save_interval = 500; % Save time by saving records only once every few iterations
save_counter = 0; % Counter to reset after saving

N_train_min = 3000; % Minimum length of training data
N_train_max = 3000; % Maximum length of training data
N_train_increment = 500; % (Minimum incr = 100) Increment value of N_train in Grid search

q_min = 20; % Min value of q in Random search
q_max = 80; % Max value of q in Random search
q_increment = 1; % Increment value of q in Grid search

p_min = 10; % Min value of p in Random search
p_max = 80; % Max value of p in Random search
p_increment = 1; % Increment value of p in Grid search

r_p_diff_min = l; % Min difference between r and p
r_p_diff_max = l; % Max difference between r and p 
r_increment = 1; % Increment value of r in Grid search     

N_train_list = N_train_min:N_train_increment:N_train_max; % List of N_train_values to search now
q_search = q_min:q_increment:q_max; % List of q parameters to search in
% p_search defined before p for loop
% r_search defined before r for loop

% Add noise once
rng('default');
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));

%% Load saved results
model_name = 'HAVOK'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
save_file = ['Data/', simulation_data_file, '_', model_name, '_error_sweep', '_sig=', sig_str, '.mat'];

try
    load(save_file);
catch
    disp('No saved results to compare to')  
    N_train_saved = [];
    MAE_saved = [];
    
end

%% Load search space records
search_space_file = ['Data/', simulation_data_file, '_', model_name, '_search_space_sig=', sig_str, '.mat'];

try
    load(search_space_file);
catch
    disp('No search space file')
    N_train_record = [];
    q_record = [];
    p_record = [];    
    search_space = zeros(0, 0, 0);
    sec_per_iteration = [0.01];
end


%% Execution time predicted
num_iterations = length(N_train_list)*length(q_search)*length(p_min:p_increment:p_max)*length(r_p_diff_min:r_increment:r_p_diff_max)
predicted_hours = mean(sec_per_iteration)*num_iterations/3600
predicted_mins = mean(sec_per_iteration)*num_iterations/60
predicted_secs = mean(sec_per_iteration)*num_iterations

%% Loop through different training lengths
for index = 1:length(N_train_list) % Loop through N_train_list
    % index is used to relate all the lists to N_train_list
    disp('--------------------------')
    N_train_timer = tic;
    N_train = N_train_list(index)
    N_train_is_new = 1; % 1 = first time using this N_train this session
    disp('searching...')

    % Get index for search space records
    N_train_i = find(N_train_record == N_train); % index of N_train in search_space
    if isempty(N_train_i) % If N_train not in record list
        N_train_i = length(N_train_record)+1; % q index becomes next in list
    end

    % Grid search for best params: p,q for each N_train
    for q = q_search
        q_is_new = 1; % 1 = first time using this N_train this session

        % Get index for search space records
        q_i = find(q_record == q); % index of q in search_space
        if isempty(q_i) % If q not in record list
            q_i = length(q_record)+1; % q index becomes next in list
        end

        w = N_train - q; % num columns of Hankel matrix
        D = (q-1)*d*Ts; % Delay duration (Dynamics in delay embedding)
        
        p_max_new = min([p_max, w, (q*m + l)]); % Max p to avoid out of bounds 
        p_search = p_min:p_increment:p_max_new; % List of p to search, for every q

        % Grid search: Search through these p values for current q
        % May not exceed hard p_max, and width of U1 and width of U2 for r
        
        for p = p_search
                       
            % Get index for search space records
            p_i = find(p_record == p); % index of p in search_space   
            if isempty(p_i) % If p not in record list
                p_i = length(p_record)+1; % q index becomes next in list
            end            
                        
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

            % Only if param combo never been searched before
            if new_search
                %% Code only executed once per N_train
                if N_train_is_new 
                    % Training data - Last sample of training is first sample of testing
                    y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
                    u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
                    t_train = t(:,end-N_test-N_train+2:end-N_test+1);

                    if ismember(N_train,N_train_saved)
                        % Use previously saved best results
                        save_index = find(N_train_saved == N_train); % Index of N_train in saved list
                        MAE_best = MAE_saved(:,save_index); % 
                        p_best = p_saved(save_index);
                        q_best = q_saved(save_index);
                        r_best = r_saved(save_index);   

                    else % If no saved data about this N_train
                        save_index = -1; % -1 means no saved data exits for N_train
                        MAE_best = Inf*ones(m,1); % MAE for every measured state
                        p_best = -1;
                        q_best = -1;
                        r_best = -1;
                    end

                    N_train_is_new = 0;
                end % if N_train_is_new

                %% Code only executed once per q
                if q_is_new 
                    timer_q = tic; % Start timer for this model evaluation

                    % Step 1: Collect and construct the snapshot matrices:
                    % According to: Discovery of Nonlinear Multiscale Systems: Sampling Strategies and Embeddings
                    % pg. 15 Delay spacing for multiscale dynamics
                    X = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
                    X2 = zeros(q*m,w); % X one step into future
                    for row = 0:q-1 % Add delay coordinates
                        X(row*m+1:(row+1)*m, :) = y_train(:, row*d + (0:w-1)*c + 1);
                        X2(row*m+1:(row+1)*m, :) = y_train(:, row*d + (0:w-1)*c + 2);
                    end

                    Upsilon = u_train(:, row*d + (0:w-1)*c + 1); % Upsilon, same indexes as last X row

                    Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

                    % Step 2: Compute the SVD of the input space Omega
                    [U1,S1,V1] = svd(Omega, 'econ');
                    %         figure, semilogy(diag(S), 'x'), hold on;
                    %         title('Singular values of Omega, showing p truncation')
                    %         plot(p,S(p,p), 'ro'), hold off;

                    % Step 3: Compute the SVD of the output space X'
                    [U2,S2,V2] = svd(X2, 'econ');
                    %         figure, semilogy(diag(S), 'x'), hold on;
                    %         title('Singular values of X2, showing r truncation')
                    %         plot(r,S(r,r), 'ro'), hold off;
                    % figure, % Plot columns of V
                    % for i=1:20    
                    %     plot(V(:,i));
                    %     pause
                    % end

                    time_q = toc(timer_q); % record time for first part in q loop

                    q_is_new = 0;
                end % if q_is_new            

                timer_p = tic; % start timer of part in p loop

                
                % Step 2.5 and 3.5: Truncate SVD matrixes with p and r
                % Do here so SVD is performed only once per q in Grid search

                % Truncate SVD matrixes of Omega
                U_tilde = U1(:, 1:p); 
                S_tilde = S1(1:p, 1:p);
                V_tilde = V1(:, 1:p);
                U1_tilde = U_tilde(1:q*m, :);
                U2_tilde = U_tilde(q*m+1:end, :);
                
                % Set r_search 
                r_search = (p - r_p_diff_max):r_increment:(p - r_p_diff_min); % List of p to search, for every q

                
                for r = r_search
                    
%                   Old line:  r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)

                    % Truncate SVD matrixes of X2
                    U_hat = U2(:, 1:r); 
                    S_hat = S2(1:r, 1:r);
                    V_hat = V2(:, 1:r);

                    % Step 4: Compute the approximation of the operators G = [A B]
                    A_tilde = U_hat'*X2*V_tilde/(S_tilde)*U1_tilde'*U_hat;
                    B_tilde = U_hat'*X2*V_tilde/(S_tilde)*U2_tilde';
        %             if (sum(abs(eig(A_tilde)) > 1) ~= 0) % If some eigenvalues are unstable due to machine tolerance
        %                  disp('Unstable eigenvalues')
        %             end

                    % If some eigenvalues are unstable due to machine tolerance,
                    % Scale them to be stable
                    count = 0;
                    while (sum(abs(eig(A_tilde)) > 1) ~= 0) && (count<10)
                        count = count+1;
                        [Ve,De] = eig(A_tilde);
                        unstable = abs(De)>1; % indexes of unstable eigenvalues
                        De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-16+count); % Normalize all unstable eigenvalues (set abs(eig) = 1)
                        A_tilde = Ve*De/(Ve); % New A with margininally stable eigenvalues
                        A_old = A_tilde;
                        A_tilde = real(A_tilde);
                    end

                    if (sum(abs(eig(A_tilde)) > 1) ~= 0) % If tilde eigenvalues are still unstable
                        % Record in search space so don't search again
                        search_space(N_train_i, q_i, p_i) = -1; 
                        N_train_record(N_train_i) = N_train; % record q in search
                        q_record(q_i) = q; % record q in search
                        p_record(p_i) = p; % record q in search
                        continue; % Exit this p loop if still unstable
                    end

                    % x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
                    % x = U_hat*x_tilde, Transform to original coordinates
                    % x_tilde = U_hat'*x, Transform to reduced order coordinates
                    % Here x is augmented state

                    A = U_hat*A_tilde*U_hat';
                    B = U_hat*B_tilde;

                    if (sum(abs(eig(A)) > 1) ~= 0) % If eigenvalues are unstable
                        % Record in search space so don't search again
                        search_space(N_train_i, q_i, p_i) = -1; 
                        N_train_record(N_train_i) = N_train; % record q in search
                        q_record(q_i) = q; % record q in search
                        p_record(p_i) = p; % record q in search
                        continue; % Exit this p loop if still unstable
                    end

                    % Time taken to train this model
                    time = time_q + toc(timer_p); % Add time taken in q loop before p chosen

                    % x_augmented(k+1) = A*x_aug(k) + B*u(k)

                    %% Compare to testing data
                    % Initial condition
                    y_hat_0 = zeros(q*m,1);
                    for row = 0:q-1 % First column of spaced Hankel matrix
                        y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)*d+1) + row*d + 1);
                    end

                    % Run model
                    Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
                    Y_hat(:,1) = y_hat_0; % Initial condition
                    for k = 1:N_test-1
                        Y_hat(:,k+1) = A*Y_hat(:,k) + B*u_test(:,k);
                    end

                    y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)

                    % Vector of Mean Absolute Error on testing data
                    MAE = sum(abs(y_hat - y_test), 2)./N_test; % For each measured state

                    % MAE metrics scaled according to range of state and average
                    % taken over all states
                    scaled_MAE = mean(MAE.*state_scale);
                    scaled_MAE_best = mean(MAE_best.*state_scale);

                    % If found best result yet, save it
                    if scaled_MAE < scaled_MAE_best
        %                 disp('Found better params:')

                        MAE_best = MAE;
                        p_best = p;
                        q_best = q;
                        r_best = r;

                        if isempty(N_train_saved) % If no saved data file existed
                            N_train_saved = [N_train];
                            MAE_saved = [MAE_best];
                            p_saved = [p_best];
                            q_saved = [q_best];
                            r_saved = [r_best];

                            save_index = 1; % save_index now has a positive value

                        elseif save_index == -1 % If first data for N_train
                            % Insert data in correct place
                            if N_train < N_train_saved(end)
                                for save_index = 1:length(N_train_saved)
                                    % Find first index where N_train is bigger
                                    if N_train < N_train_saved(save_index)
                                        % Insert at current saved_index
                                        N_train_saved = insert(N_train_saved, save_index, N_train);
                                        MAE_saved = insert(MAE_saved, save_index, MAE_best);
                                        p_saved = insert(p_saved, save_index, p_best);
                                        q_saved = insert(q_saved, save_index, q_best);
                                        r_saved = insert(r_saved, save_index, r_best);            
                                        break; 
                                    end
                                end

                            else % If N_train is the biggest yet, add to end of list
                                N_train_saved = [N_train_saved, N_train];
                                MAE_saved = [MAE_saved, MAE_best];
                                p_saved = [p_saved, p_best];
                                q_saved = [q_saved, q_best];
                                r_saved = [r_saved, r_best];

                                save_index = length(N_train_saved); % save_index now has a positive value
                            end


                        else % Replace previous saved data for N_train

                            MAE_saved(:,save_index) = MAE_best;
                            p_saved(save_index) = p_best;
                            q_saved(save_index) = q_best;
                            r_saved(save_index) = r_best;

                        end

                        % Save results
                        save(save_file, 'N_train_saved', 'MAE_saved', 'p_saved', 'q_saved', 'r_saved');
                    end

                    %% Record in search space 
                    % (Also happens when break from unstable eigenvalues)
                    search_space(N_train_i, q_i, p_i) = scaled_MAE;
                    N_train_record(N_train_i) = N_train; % record q in search
                    q_record(q_i) = q; % record q in search
                    p_record(p_i) = p; % record q in search

                    % Save every few iterations to keep progress, yet save time
    %                 save(search_space_file, 'search_space', 'N_train_record', 'q_record', 'p_record', 'sec_per_iteration')
                    if save_counter > save_interval
                        disp('save')
                        save(search_space_file, 'search_space', 'N_train_record', 'q_record', 'p_record', 'sec_per_iteration')
                        save_counter = 0;
                    end
                    save_counter = save_counter+1;
                
                end % p loop  
                
            else
                % Do nothing, param combo searced before
%                 disp('Params were searched before')
            end  % if(new_search)
            

        end % p loop
    
    end % q loop
   
%     disp('Results:')
%     MAE_best
%     p_best
%     q_best
%     time
    N_train_time = toc(N_train_timer);
    disp(['Spent ', num2str(N_train_time), ' seconds searching this N_train'])
    disp('--------------------------')

end

%% Save search_space 
total_time = toc(total_timer);
sec_per_iteration = [sec_per_iteration; total_time/num_iterations];
save(search_space_file, 'search_space', 'N_train_record', 'q_record', 'p_record', 'sec_per_iteration')

%% Plot saved results
disp('Plotting...')
close all;

for fig_num = 1:m % Plot error of every measured state
    figure(fig_num), semilogy(N_train_saved,MAE_saved(fig_num,:),'x'), title(['MAE of x(', num2str(fig_num), ') vs N-train'])
end
fig_num = fig_num + 1; % Increment fig_num for next figure
figure(fig_num), plot(N_train_saved,p_saved,'x'), title('p vs N-train')
fig_num = fig_num + 1; 
figure(fig_num), plot(N_train_saved,q_saved,'x'), title('q vs N-train')
fig_num = fig_num + 1;
figure(fig_num), plot(N_train_saved,r_saved,'x'), title('time vs N-train')
fig_num = fig_num + 1;

figure(fig_num)
fig_num = fig_num + 1;
p_table = tabulate(p_saved);
subplot(1,2,1), bar(p_table(:,1), p_table(:,2)), title('Frequency of p');
q_table = tabulate(q_saved);
subplot(1,2,2), bar(q_table(:,1), q_table(:,2)), title('Frequency of q');


toc(total_timer);
disp('--------------------------')
disp('END of HAVOK N_train sweep')

%% Local functions
function X_hat = plot_model(A,B,U_data,t,x0)
    N = max(size(t));   % Number of time steps 
    X_hat = zeros(length(x0),N); % Estimated X from model
    X_hat(:,1) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
    plot(t, X_hat);     % Estimated x
end

function X_hat = run_model(A,B,U_data,t,x0)
    N = max(size(t));   % Number of time steps 
    X_hat = zeros(length(x0),N); % Estimated X from model
    X_hat(:,1:size(x0)) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
end

function new_array = insert(array, index, entry)
    if index == 1 % to avoid index-1 = 0
        new_array = [entry, array];
    else
        new_array = [array(:, 1:index-1), entry, array(:, index:end)];
    end
end






