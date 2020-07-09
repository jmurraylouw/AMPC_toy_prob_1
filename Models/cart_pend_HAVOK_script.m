%% HAVOK with control - of cart pendulum
% Estimate linear model from data
% Partial state feedback

%% Try to save previous results of random search and continue with them

%% Read data
tic;

close all;
load('cartpend_random_1.mat') % Load simulation data
% x0 = [1; -0.2; -0.5; 0.8]
u_data  = out.u.Data';
x_data  = out.x.Data';
y_data  = x_data([1,3],:); % Measurement data (x and theta)

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
N_test = 5000; % Num of data samples for testing
c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)
% N_train; % Num of data samples for training, rest for testing
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
N_train_min = 10/Ts; % Minimum length of training data
N_train_max = 50/Ts; % Maximum length of training data
max_iterations = 1000; % Maximum number of iterations allowed in Random search
% p = 34; % Truncated rank of system
p_min = 20; % Min value of p in Random search
p_max = 40; % Max value of p in Random search
% q = 1000; % number of delays
q_min = 50; % Min value of q in Random search
% q_max; % Max value of q in Random search

% Add noise once
y_data_noise = y_data + sigma*randn(size(y_data));

% Lists to save results for different N_train
N_train_list = N_train_min:1000:N_train_max; % List of N_train to compute model for
MAE_list = Inf + zeros(m,length(N_train_list)); % Save best MAE for each N_train
p_list = -1 + zeros(1,length(N_train_list)); % Save p of best MAE for each N_train
q_list = -1 + zeros(1,length(N_train_list)); % Save q of best MAE for each N_train
time_list = -1 + zeros(1,length(N_train_list)); % Save time taken for each N_train

% Loop through different training lengths
for index = 1:length(N_train_list) % Loop through N_train_list
    % index is used to relate all the lists to N_train_list
    N_train = N_train_list(index); 
    
    p_best = NaN;
    q_best = NaN;
    MAE_best = Inf*[1;1];

%     q_max = floor(N_train/2); % Max q when Hankel is a square
%     for q = 
    Random search for best hyperparameters
    for iteration = 1:max_iterations % Loop truncateed rank
        timer = tic; % Start timer for this model evaluation
        
        q_max = floor(N_train/2); % Max q when Hankel is a square
        q = randi([q_min, q_max]); % Scaled to get better uniform distribution
        p = randi([p_min, p_max]);
        
        w = N_train - q; % num columns of Hankel matrix
        
        % numel_H = q*m*w;
        % log10(numel_H)
        % time_predict = 6e-6*numel_H

        D = (q-1)*d*Ts; % Delay duration

        % Testing data - Last 50 s is for testing and one sample overlaps training 
        y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
        u_test = u_data(:,end-N_test+1:end);
        t_test = t(:,end-N_test+1:end);

        % Training data - Last sample of training is first sample of testing
        y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
        u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
        t_train = t(:,end-N_test-N_train+2:end-N_test+1);

        r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)

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

        % Step 2: Compute and truncate the SVD of the input space Omega
        [U,S,V] = svd(Omega, 'econ');
%         figure, semilogy(diag(S), 'x'), hold on;
%         title('Singular values of Omega, showing p truncation')
%         plot(p,S(p,p), 'ro'), hold off;
        % figure, % Plot columns of V
        % for i=1:20    
        %     plot(V(:,i));
        %     pause
        % end
        U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
        S_tilde = S(1:p, 1:p);
        V_tilde = V(:, 1:p);
        U1_tilde = U_tilde(1:q*m, :);
        U2_tilde = U_tilde(q*m+1:end, :);

        % Step 3: Compute the SVD of the output space X'
        [U,S,V] = svd(X2, 'econ');
%         figure, semilogy(diag(S), 'x'), hold on;
%         title('Singular values of X2, showing r truncation')
%         plot(r,S(r,r), 'ro'), hold off;
        U_hat = U(:, 1:r); % Truncate SVD matrixes of X2
        S_hat = S(1:r, 1:r);
        V_hat = V(:, 1:r);

        % Step 4: Compute the approximation of the operators G = [A B]
        A_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U1_tilde'*U_hat;
        B_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U2_tilde';
        if (sum(abs(eig(A_tilde)) > 1) ~= 0) % If some eigenvalues are unstable due to machine tolerance
%             disp('Unstable eigenvalues')
        end

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
        assert(sum(abs(eig(A_tilde)) > 1) == 0, 'Unstable eigenvalues'); % Check if all eigenvalues are stable (magnitude <= 1)

        % x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
        % x = U_hat*x_tilde, Transform to original coordinates
        % x_tilde = U_hat'*x, Transform to reduced order coordinates
        % Here x is augmented state

        A = U_hat*A_tilde*U_hat';
        B = U_hat*B_tilde;

        if (sum(abs(eig(A)) > 1) ~= 0) % If some eigenvalues are unstable due to machine tolerance
            disp('Still unstable eigenvalues')
            break;
        end

        time = toc(timer);
        
        % x_augmented(k+1) = A*x_aug(k) + B*u(k)
        % Ignore eigenmodes Step 5 and 6

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

        if mean(MAE) < mean(MAE_best)
            MAE_best = MAE;
            p_best = p;
            q_best = q;
        end

        % %% Compare to training data
        % disp(6)
        % % Initial conditions
        % y_hat_02 = zeros(q*m,1);
        % for row = 0:q-1 % Create first column of spaced Hankel matrix
        %     y_hat_02(row*m+1:(row+1)*m, 1) = y_train(:, row*d + 1);
        % end
        % k_start = row*d + 1; % First k to start at
        % 
        % Y_hat2 = zeros(length(y_hat_0),N_train); % ??? Estimated X from model
        % Y_hat2(:,k_start) = y_hat_02; % Initial conditions, insert at first k
        % for k = k_start:N_train-1
        %     Y_hat2(:,k+1) = A*Y_hat2(:,k) + B*u_train(:,k);
        % end
        % y_hat2 = Y_hat2(end-m+1:end, :); % Extract only non-delay time series (last m rows)
        % 
        % disp('Run model on training data')

        % %% Plot data vs model
        % figure;
        % plot(t, y_data); 
        % hold on;
        % plot(t, u_data, ':', 'LineWidth', 1);
        % plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
        % plot(t_train, y_hat2, '--', 'LineWidth', 1); % Plot only non-delay coordinate  
        % plot((D + t(N-N_test-N_train)).*[1,1], ylim, 'r');
        % plot(t(N-N_test-N_train).*[1,1], ylim, 'r');
        % plot(t(N-N_test).*[1,1], ylim, 'k');
        % title('Training and Testing data vs Model');
        % % legend('x', 'theta', 'input', 'x_hat', 'theta_hat', 'D', 't(final sample)')
        % hold off;
        % 
        % toc;

    end

    MAE_list(:,index) = MAE_best; % Save best MAE related to current N_train
    p_list(:,index) = p_best;
    q_list(:,index) = q_best;
    time_list(:,index) = time; % Save time taken to compute model for each N_train
      
    disp('--------------------------')
    N_train
    MAE_best
    p_best
    q_best
    time
    
end

%% Merge new and saved results
% If two results are for the same N_train, keep the one with the smallest error
try
    load('Data\N_train_error_time_HAVOK_sig=0.mat');
catch
    disp('No saved results to compare to')
    N_train_saved = N_train_list;
    MAE_saved = MAE_list;
    p_saved = p_list;
    q_saved = q_list;
    time_saved = time_list;
end
% 
% saved = [N_train_saved', mean(MAE_saved',2)]
% list = [N_train_list', mean(MAE_list', 2)]

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

%%
figure(1), plot(N_train_saved,MAE_saved(1,:)), title('MAE x vs N-train')
figure(2), plot(N_train_saved,MAE_saved(2,:)), title('MAE theta vs N-train')
figure(3), plot(N_train_saved,p_saved), title('p vs N-train')
figure(4), plot(N_train_saved,q_saved), title('q vs N-train')
figure(5), plot(N_train_saved,time_saved), title('time vs N-train')

toc;
disp('------------')
disp('END of HAVOK')


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






