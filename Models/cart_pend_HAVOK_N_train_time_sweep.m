%% HAVOK N_train vs time sweep
% HAVOK with control - of cart pendulum
% Calculate average time taken to generate model with different training data sizes

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

time_iterations = 100; % Num iterations per model to average time
sigma = 0.0; % Noise standard deviation
c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)
% N_train; % Num of data samples for training, rest for testing
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
% p = Truncated rank of system of Omega
% r = Truncated rank of system of X2
% q = number of delays

% Add noise once
rng('default');
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));


%% Load saved results

model_name = 'HAVOK'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
load_file = ['Data\', model_name, '_N_train_vs_error', '_sig=', sig_str, '.mat'];
save_file = ['Data\', model_name, '_N_train_vs_time', '_sig=', sig_str, '.mat'];
%%
try
    load(load_file);   
catch
    error('No saved results to load')  
end

% Number of iterations predicted
predicted_iterations = time_iterations*length(N_train_saved)

%% Loop through different training lengths
% Empty list to store average time per N_train
time_list = zeros(1,length(N_train_saved));
for index = 1:length(N_train_saved)
    N_train = N_train_saved(index)
    p = p_saved(index);
    q = q_saved(index);
            
    r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)
    w = N_train - q; % num columns of Hankel matrix
    D = (q-1)*d*Ts; % Delay duration (Dynamics in delay embedding)

    % Training data - Last sample of training is first sample of testing
    y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
    u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
    t_train = t(:,end-N_test-N_train+2:end-N_test+1);

    % Empty time list to average for this N_train
    time_list_i = zeros(1,length(time_iterations));
    for i = 1:time_iterations

        timer_model = tic; % Start timer for this model evaluation

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

        % Step 3: Compute the SVD of the output space X'
        [U2,S2,V2] = svd(X2, 'econ');

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
            error('eigenvalues are unstable'); % Exit this p loop if still unstable
        end

        A = U_hat*A_tilde*U_hat';
        B = U_hat*B_tilde;

        if (sum(abs(eig(A)) > 1) ~= 0) % If eigenvalues are unstable
            error('eigenvalues are unstable');
        end

        % Time taken to train this model
        time_list_i(i) = toc(timer_model); % Add time taken in q loop before p chosen

    end % End of iterations for this N_train
    
%     mean_time = mean(time_list_i)
%     disp('--------------------')
%     plot(time_list_i)
%     title(num2str(N_train))

    time_list(index) = mean(time_list_i);
end % End of N_train loop

plot(time_list)
toc(total_timer);

save(save_file, 'time_list');

disp('-----------------------')
disp('END of HAVOK time sweep')


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






