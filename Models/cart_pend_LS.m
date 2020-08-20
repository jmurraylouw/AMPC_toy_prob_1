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

% Training data - Last sample of training is first sample of testing
N_train = 2000; % Num of data samples for training, rest for testing
x_train = x_data(:,end-N_test-N_train+2:end-N_test+1);
y_train = y_data(:,end-N_test-N_train+2:end-N_test+1);
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

% Add noise
sigma = 0.1;
y_train = y_train + sigma*randn(size(y_train));

% Best delays for N_train = 2000 (found with cartpend_LS_delay_error_sweep.m)
switch sigma
    case 0.1
        delays = 85;
    case 0.01
        delays = 66;
    case 0
        delays = 2;
end

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

%% Run model on testing data
t_hat = t_test(:, delays:end-1); % Shorten to fit with x_hat, due to initial conditions
u_hat = u_test(:, delays:end-1);
y_hat = zeros(na, length(t_hat));

y_hat_0 = []; % Initial condition
for i = 1:delays
    y_hat_0 = [y_test(:,i); y_hat_0];
end
y_hat(:,1) = y_hat_0; % Initial condition

for index = 1:length(t_hat)-1
    y_hat(:,index+1) = A*y_hat(:,index) + B*u_hat(index);
end

y_hat = y_hat(1:ny, :); % Only keep non-delay coordinates

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test(:, delays:end-1)), 2)./size(y_hat,2); % For each measured state
MAE

%%
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





