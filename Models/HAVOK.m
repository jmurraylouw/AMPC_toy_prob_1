% Implentation of Hankel Alternative View Of Koopman

close all
rng('default');
rng(1); % Repeatable random numbers

% Simulation parameters
Ts = 0.01; % Sample time
t = 0:Ts:100; % Simulation time period
N = length(t);
x0 = [0.5; 0]; % Initial conditions
n = length(x0);
f = @nl_msd; % Model for this simulation

% Input data
u_data = zeros(1, length(t));
rand_u = 0;
rand_interval = 2; % Size of
sigma_u = 1;

for i = 1:length(t)
    if mod(i, rand_interval/Ts) == 0 % Only take new random input value every few seconds
        rand_u = sigma_u*randn(1,1);
    end
    u_data(:,i) = rand_u;
end

% Run simulation
% [t,x] = ode45(@(t,x) nl_msd(x,0), tspan, x0); % simulate with no input

% Solve for small intervals with constant u
x = zeros(n,N); % Empty prediction matrix
x(:,1) = x0; 

for i = 1:N - 1
    x00 = x(:,i); % Initial condition for this timestep
    u = u_data(:,i); %(U_test(i,:) + U_test(i+1,:))/2; % Assume constant u at average of time interval
    [t_1,x_1] = ode45(@(t_1,x_1) f(x_1,u), t(i:i+1), x00);
    x(:,i+1) = x_1(end,:)';
end

% Extract data
u_data  = u_data;
x_data  = x';
y_data  = x_data([1,2],:); % Measurement data (x, z, theta)
t       = t';

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 5000; % Num of data samples for testing
x_test = x_data(:, (end-N_test+1):end);
y_test = y_data(:, (end-N_test+1):end); % One sample of testing data overlaps for initial condition
u_test = u_data(:, (end-N_test+1):end);
t_test =      t(:, (end-N_test+1):end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% Add noise
sigma = 0; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data - Last sample of training is first sample of testing
N_train = 5000; % Number of sampels in training data
y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

% Parameters
q = 400;
p = 300;
w = N_train - q + 1; % num columns of Hankel matrix
D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

% Create Hankel matrix with measurements
Y = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
for row = 0:q-1 % Add delay coordinates
    Y(row*m+1:(row+1)*m, :) = y_train(:, row + (0:w-1) + 1);
end

% SVD of the Hankel matrix
[U1,S1,V1] = svd(Y, 'econ');
% figure, semilogy(diag(S1), 'x'), hold on;
% title('Singular values of Omega, showing p truncation')
% plot(p,S1(p,p), 'ro'), hold off;

% Truncate SVD matrixes
U_tilde = U1(:, 1:p); 
S_tilde = S1(1:p, 1:p);
V_tilde = V1(:, 1:p);

% Setup V2 one timestep into future from V1
V_til_2 = V_tilde(2:end  , :)'; % Turnd on side (wide short matrix)
V_til_1 = V_tilde(1:end-1, :)';

% Based on DMD control example video by Steve Brunton
U = u_train(:, q:end-1); % Leave out last time step to match V_til_1
VU = [V_til_1; U]; % Combined matrix of V and U, above and below
AB = V_til_2*pinv(VU); % combined A and B matrix, side by side
A_tilde  = AB(:,1:p); % Extract A matrix
B_tilde  = AB(:,(p+1):end);

% DMD on V
% A_tilde = V_til_2*pinv(V_til_1); % Matrix to propogate V' forward in time. Note transpose to turn V into fat/horizontal matrix

% convert to x coordinates
A = (U_tilde*S_tilde)*A_tilde*pinv(U_tilde*S_tilde);

% DMD
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

A = Y2*pinv(Y1);

%% Compare to testing data

%% Run with A_tilde and v
figure;
plot(V1(:,1:5))

% Initial condition
v_hat_0 = V_tilde(end,:)';

% Run model
V_hat = zeros(length(v_hat_0),N_test); % Empty estimated Y
V_hat(:,1) = v_hat_0; % Initial condition
for k = 1:N_test-1
    V_hat(:,k+1) = A_tilde*V_hat(:,k) + B_tilde*u_test(:,k);
end

Y_hat2 = (U_tilde*S_tilde)*V_hat; % Convert to Y
y_hat2 = Y_hat2(end-m+1:end, :); % Extract only non-delay time series (last m rows)

% Vector of Mean Absolute Error on testing data
MAE_til = sum(abs(y_hat2 - y_test), 2)./N_test % For each measured state


%% Run with A and x

% Initial condition
y_hat_0 = zeros(q*m,1);
for row = 0:q-1 % First column of spaced Hankel matrix
    y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)+1) + row + 1);
end

% Run model
Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
Y_hat(:,1) = y_hat_0; % Initial condition
for k = 1:N_test-1
    Y_hat(:,k+1) = A*Y_hat(:,k);
end

y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

%% Compare MAE and MAE_til
MAE_error_percent = (MAE - MAE_til)./MAE_til*100

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
plot((D + t(N-N_test-N_train)).*[1,1], ylim, 'r');
plot(t(N-N_test-N_train).*[1,1], ylim, 'k');
plot(t(N-N_test).*[1,1], ylim, 'k');
title('Training and Testing data vs Model');
hold off;


function dx = nl_msd(x,u)
    % Non-linear mass spring damper
    dx = zeros(2,1);

    % Parameters
    m = 1;
    k = 1.1;
    k_nl = 0.8;
    c = 0.05;
    
    % State space ODE
    dx(1,1) = 1/m*(x(2));
    dx(2,1) = 1/m*(-c*x(2) - k*x(1) - k_nl*abs(x(1)) + u);
end
