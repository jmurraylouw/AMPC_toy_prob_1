%% OLS of cart pendulum
% Estimate System matrixes with moving window of data in real time
% Partial state feedback

% Setup
rng(0);
random_generator = rng; % repeatable random numbers

% Read simulation data
load('cartpend_random_1');
x_data = out.x.Data'; % each row is timeseries of a state
y_data = x_data([1,3],:);
u_data = out.u.Data';
t = out.x.Time';
Ts = t(2)-t(1);
N = length(t); % Total number of samples in data
nx = size(x_data,1);
ny = size(y_data,1);
nu = size(u_data,1);

% Testing data - Last 50 s is for testing and one sample overlaps training
N_test = 5000; % Num of data samples for testing
% One sample of testing data overlaps for initial condition
x_test = x_data(:,end-N_test+1:end); % each row is timeseries of a state
y_test = y_data(:,end-N_test+1:end); 
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Add noise
sigma = 0.01;
y_data_noise = y_data + sigma*randn(size(y_data));

% Best delays for N_train = 2000 (found with cartpend_LS_delay_error_sweep.m)
switch sigma
    case 0.1
        delays = 85;
    case 0.01
        delays = 66;
    case 0
        delays = 2;
end

% List of N_train values to sweep through
N_train_min = 100; % Minimum length of training data
N_train_max = 5000; % Maximum length of training data
N_train_increment = 100; % (Minimum incr = 100) Increment value of N_train in Grid search

N_train_saved = N_train_min:N_train_increment:N_train_max; % List of N_train_values to search now

MAE_saved = zeros(ny,length(N_train_saved))-1; % Save MAE for corresponding N_train value
time_saved  = zeros(1,length(N_train_saved))-1; % Save time taken to generate model for N_train values

for index = 1:length(N_train_saved) % Loop through N_train_list
    % index is used to relate all the lists to N_train_list
    N_train = N_train_saved(index) % Current N_train value to evaluate
    
    % Training data - Last sample of training is first sample of testing
    y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1);
    u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
    t_train = t(:,end-N_test-N_train+2:end-N_test+1);

    timer = tic;
    % Create delays matrix
    % Note that X is made from measurements, not full state
    X = []; % Augmented state with delay coordinates [Y(k); Y(k-1); Y(k-2); ...]
    for i = 1:delays
        X = [y_train(:, i:end-delays+i); X];  
    end

    na = size(X,1); % Length of augmented state vector

    X2 = X(:, 2:end); % X advanced 1 step into future
    X = X(:, 1:end-1); % Cut off last sample

    Upsilon = u_train(:, delays:end-1); % U snapshot

    Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

    % Least squares
    AB = X2*pinv(Omega);

    A = AB(:,1:na);
    B = AB(:,na+1:end);

    % Record time
    time_saved(index) = toc(timer); % Save time taken
    
    %% Run model on testing data
    t_hat = t_test(:, delays:end-1); % Shorten to fit with x_hat, due to initial conditions
    u_hat = u_test(:, delays:end-1);
    y_hat = zeros(na, length(t_hat));

    y_hat_0 = []; % Initial condition
    for i = 1:delays
        y_hat_0 = [y_test(:,i); y_hat_0];
    end
    y_hat(:,1) = y_hat_0; % Initial condition

    for i = 1:length(t_hat)-1
        y_hat(:,i+1) = A*y_hat(:,i) + B*u_hat(i);
    end

    y_hat = y_hat(1:ny, :); % Only keep non-delay coordinates

    % Vector of Mean Absolute Error on testing data
    MAE = sum(abs(y_hat - y_test(:, delays:end-1)), 2)./size(y_hat,2); % For each measured state
    
    % Save MAE value for current N_train
    MAE_saved(:,index) = MAE;
end

%% Save results
model_name = 'LS'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
save_file = ['Data\', model_name, '_N_train_vs_error', '_sig=', sig_str, '.mat'];

save(save_file, 'N_train_saved', 'MAE_saved', 'delays', 'time_saved');

%% Load results (if did not run with this sigma again)
load(save_file)

%% Plot time taken vs N_train
figure()
semilogy(N_train_saved, time_saved)
title('N-train vs time')

%% Plot N-train vs MAE
figure()
semilogy(N_train_saved, MAE_saved)
title('N-train vs MAE')

%% plot testing data
figure()
hold on;
title('Testing data')
plot(t_test, y_test)
plot(t_hat, y_hat, '--')
stop

%% Run model on training data
t_hat = t_train(:, delays:end-1); % Shorten to fit with x_hat, due to initial conditions
u_hat = u_train(:, delays:end-1);
y_hat = zeros(na,length(t_hat));

y_hat_0 = []; % Initial condition
for i = 1:delays
    y_hat_0 = [y_train(:,i); y_hat_0];
end
y_hat(:,1) = y_hat_0; % Initial condition

for index = 1:length(t_hat)-1
    y_hat(:,index+1) = A*y_hat(:,index) + B*u_hat(index);
end

y_hat = y_hat(1:ny, :); % Only keep non-delay coordinates

% Vector of Mean Absolute Error on testing data
MAE_train = sum(abs(y_hat - y_train(:, delays:end-1)), 2)./size(y_hat,2) % For each measured state

figure()
hold on;
title('Training data')
plot(t_train, y_train)
plot(t_hat, y_hat, '--')

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

function X_hat = run_model(A,B,U_data,x0)
    N = max(size(U_data));   % Number of time steps 
    X_hat = zeros(length(x0),N); % Estimated X from model
    X_hat(:,1) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
end





