%% HAVOK with control - of cart pendulum
% Estimate linear model from data
% Partial state feedback

%% Try to save previous results of random search and continue with them

%% Read data
total_timer = tic;

close all;
load('cartpend_random_1.mat') % Load simulation data
% x0 = [1; -0.2; -0.5; 0.8]
u_data  = out.u.Data';
x_data  = out.x.Data';
y_data  = x_data([1,3],:); % Measurement data (x and theta)
t       = out.tout'; % Time

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 5000; % Num of data samples for testing
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
t  = out.tout';
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% [X_p,Y_delays] = meshgrid(1:delays(end), 1:delays(end)); % Values for surface plot
% RMSE_matrix = zeros(delays(end), delays(end)); % Empty matrix of errors

%% Parameters
% Very dependant on choice of p, r, q

sigma = 0.0; % Noise standard deviation
c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)
% N_train; % Num of data samples for training, rest for testing
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
N_train_min = 1/Ts; % Minimum length of training data
N_train_max = 50/Ts; % Maximum length of training data
% p = Truncated rank of system of Omega
% r = Truncated rank of system of X2
% q = number of delays
p_min = 8; % Min value of p in Random search
p_max = 8; % Max value of p in Random search
q_min = 4; % Min value of q in Random search
q_max = 4; % Max value of q in Random search
q_increment = 1; % Increment value of q in Grid search
p_increment = 1; % Increment value of p in Grid search
N_train_increment = 10; % Increment value of N_train in Grid search

N_train_list = N_train_min:N_train_increment:N_train_max;
q_search = q_min:q_increment:q_max; % List of q parameters to search in

% Add noise once
rng('default');
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));

% Number of iterations predicted
num_iterations = (p_max - p_min)/p_increment*((N_train_max+N_train_min)/4 - q_min)/q_increment*length(N_train_list)

%% Load saved results

model_name = 'HAVOK'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
save_file = ['Data\N_train_error_time_', model_name, '_sig=', sig_str, '.mat'];

try
    load(save_file);
catch
    disp('No saved results to compare to')  
    N_train_saved = [];
end

%% Loop through different training lengths
for index = 1:length(N_train_list) % Loop through N_train_list
    % index is used to relate all the lists to N_train_list
    disp('--------------------------')
    N_train = N_train_list(index)
    disp('searching...')

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
        time_best = time_saved(save_index);   

    else % If no saved data about this N_train
        save_index = -1; % -1 means no saved data exits for N_train
        MAE_best = Inf*[1;1];
        p_best = NaN;
        q_best = NaN;
        time_best = NaN;
    end

%     q_max = floor(N_train/3); % Max q when Hankel is a square
    
    % Grid search for best params: p,q for each N_train
    for q = q_search

%     %  Random search for best hyperparameters
%     for iteration = 1:max_iterations % Loop truncateed rank         
%         q_max = floor(N_train/2); % Max q when Hankel is a square
%         q = randi([q_min, q_max]); % Scaled to get better uniform distribution
%         p = randi([p_min, p_max]);

        timer_q = tic; % Start timer for this model evaluation
        w = N_train - q; % num columns of Hankel matrix

        D = (q-1)*d*Ts; % Delay duration (Dynamics in delay embedding)

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
        
        % Grid search: Search through these p values for current q    
        p_search = p_min:p_increment:min([p_max,q*m]); % List of p to search, for every q
        for p = p_search
            
            timer_p = tic; % start timer of part in p loop
            
            r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)

            % Step 2.5 and 3.5: Truncate SVD matrixes with p and r
            % Do here so SVD is performed only once per q in Grid search
            
            % Truncate SVD matrixes of Omega
            U_tilde = U1(:, 1:p); 
            S_tilde = S1(1:p, 1:p);
            V_tilde = V1(:, 1:p);
            U1_tilde = U_tilde(1:q*m, :);
            U2_tilde = U_tilde(q*m+1:end, :);

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
            while (sum(abs(eig(A_tilde)) > 1) ~= 0) 
                count = count+1;
                [Ve,De] = eig(A_tilde);
                unstable = abs(De)>1; % indexes of unstable eigenvalues
                De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-16+count); % Normalize all unstable eigenvalues (set abs(eig) = 1)
                A_tilde = Ve*De/(Ve); % New A with margininally stable eigenvalues
                A_old = A_tilde;
                A_tilde = real(A_tilde);
                if(count>10)
                    'break'
                    break
                end
            end
            
            if (sum(abs(eig(A_tilde)) > 1) ~= 0) % If tilde eigenvalues are still unstable
                break; % Exit this p loop if still unstable
            end
            
            % x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
            % x = U_hat*x_tilde, Transform to original coordinates
            % x_tilde = U_hat'*x, Transform to reduced order coordinates
            % Here x is augmented state

            A = U_hat*A_tilde*U_hat';
            B = U_hat*B_tilde;

            if (sum(abs(eig(A)) > 1) ~= 0) % If eigenvalues are unstable
                break; % Exit this p loop if still unstable
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

            % If found best result yet, save it
            if mean(MAE) < mean(MAE_best)
                disp('Found better params:')
                MAE_best = MAE;
                p_best = p
                q_best = q
                time_best = time;
                
                if isempty(N_train_saved) % If no saved data file existed
                    N_train_saved = [N_train];
                    MAE_saved = [MAE_best];
                    p_saved = [p_best];
                    q_saved = [q_best];
                    time_saved = [time_best];
                    
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
                                time_saved = insert(time_saved, save_index, time_best);            
                                break; 
                            end
                        end
                        
                    else % If N_train is the biggest yet, add to end of list
                        N_train_saved = [N_train_saved, N_train];
                        MAE_saved = [MAE_saved, MAE_best];
                        p_saved = [p_saved, p_best];
                        q_saved = [q_saved, q_best];
                        time_saved = [time_saved, time_best];
                        
                        save_index = length(N_train_saved); % save_index now has a positive value
                    end

                    
                else % Replace previous saved data for N_train
                
                    MAE_saved(:,save_index) = MAE_best;
                    p_saved(save_index) = p_best;
                    q_saved(save_index) = q_best;
                    time_saved(save_index) = time_best;
                    
                end
                
                % Save results
                save(save_file, 'N_train_saved', 'MAE_saved', 'p_saved', 'q_saved', 'time_saved');
            end

        end % end of p loop
    
    end % end of q loop
   
    disp('Results:')
    MAE_best
    p_best
    q_best
    time
    disp('--------------------------')
    
end


%% Plot saved results
disp('Plotting...')
close all;
figure(1), semilogy(N_train_saved,MAE_saved(1,:)), title('MAE x vs N-train')
figure(2), semilogy(N_train_saved,MAE_saved(2,:)), title('MAE theta vs N-train')
figure(3), plot(N_train_saved,p_saved), title('p vs N-train')
figure(4), plot(N_train_saved,q_saved), title('q vs N-train')
figure(5), plot(N_train_saved,time_saved), title('time vs N-train')

figure(6)
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






