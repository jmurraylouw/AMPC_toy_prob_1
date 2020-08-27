%% HAVOK with control - of cart pendulum
% Estimate linear model from data
% Partial state feedback

%% Try to save previous results of random search and continue with them

%% Read data

close all;
% clear all;

total_timer = tic;

% load('cartpend_random_1.mat') % Load simulation data
load('Data/cartpend_disturbance_and_PID_1.mat') % Load simulation data

% Extract data
u_data  = out.u.Data';
x_data  = out.x.Data';
y_data  = x_data([1,3],:); % Measurement data (x, z, theta)
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
u_bar = [-5]; % Mean input needed to keep at a fized point
% u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 5000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end)
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

sigma = 0.01; % Noise standard deviation
N_train = 4000; % Number of sampels in training data
c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
% p % Truncated rank of system of Omega
% q % number of delays
% r = Truncated rank of system of X2

% Add noise
rng('default');
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));

%% Load saved results
model_name = 'HAVOK'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string

save_file = ['Data/', model_name, '_N_train_vs_error', '_sig=', sig_str, '.mat'];

try
    load(save_file);
    if ismember(N_train,N_train_saved)
        % Use previously saved best results
        save_index = find(N_train_saved == N_train); % Index of N_train in saved list
        MAE = MAE_saved(:,save_index) 
        p = p_saved(save_index)
        q = q_saved(save_index)
        time = time_saved(save_index)
        
%         % Override
%         disp('Override')
%         disp('------------------')
% 
%         p = 15

    else
        N_train
        error('No saved results for this N_train value')
        
%         disp('Override')
%         q = 4
%         p = 8
    end
    
catch
    disp('Saved results file does not exist')  
    q = 50
    p = 30
    
end
        
r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)
w = N_train - q; % num columns of Hankel matrix
D = (q-1)*d*Ts; % Delay duration (Dynamics in delay embedding)

% Training data - Last sample of training is first sample of testing
y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

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

Upsilon = u_train(:, row*d + (0:w-1)*c + 1) - u_bar; % Upsilon, same indexes as last X row

Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

% Step 2: Compute the SVD of the input space Omega
[U1,S1,V1] = svd(Omega, 'econ');
figure, semilogy(diag(S1), 'x'), hold on;
title('Singular values of Omega, showing p truncation')
plot(p,S1(p,p), 'ro'), hold off;

% Step 3: Compute the SVD of the output space X'
[U2,S2,V2] = svd(X2, 'econ');
figure, semilogy(diag(S2), 'x'), hold on;
title('Singular values of X2, showing r truncation')
plot(r,S2(r,r), 'ro'), hold off;
% figure, % Plot columns of V
% for i=1:20    
%     plot(V(:,i));
%     pause
% end

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
    error('eigenvalues are unstable'); % Exit this p loop if still unstable
end

% x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
% x = U_hat*x_tilde, Transform to original coordinates
% x_tilde = U_hat'*x, Transform to reduced order coordinates
% Here x is augmented state

A = U_hat*A_tilde*U_hat';
B = U_hat*B_tilde;

if (sum(abs(eig(A)) > 1) ~= 0) % If eigenvalues are unstable
    error('eigenvalues are unstable');
end

% Time taken to train this model
time = toc(timer_model) % Add time taken in q loop before p chosen

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
    Y_hat(:,k+1) = A*Y_hat(:,k) + B*(u_test(:,k) - u_bar);
end

y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state


%% Compare to training data
disp(6)
% Initial conditions
y_hat_02 = zeros(q*m,1);
for row = 0:q-1 % Create first column of spaced Hankel matrix
    y_hat_02(row*m+1:(row+1)*m, 1) = y_train(:, row*d + 1);
end
k_start = row*d + 1; % First k to start at

Y_hat2 = zeros(length(y_hat_0),N_train); % ??? Estimated X from model
Y_hat2(:,k_start) = y_hat_02; % Initial conditions, insert at first k
for k = k_start:N_train-1
    Y_hat2(:,k+1) = A*Y_hat2(:,k) + B*(u_train(:,k) - u_bar);
end
y_hat2 = Y_hat2(end-m+1:end, :); % Extract only non-delay time series (last m rows)

disp('Run model on training data')

%% Compare to linearised model

f = @cartpend;
[A_lin, B_lin] = linearise_floating_pend_2D(f,zeros(n,1),u_bar) % Get linearised model

% Initial condition
x_hat3 = zeros(size(x_test)); % Empty estimated y from pre-determined linearised model
x_hat3(:,1) = x_test(:,1); % Initial condition

% Run model
% Solve for small intervals with constant u
for i=1:N_test-1
    x00 = x_hat3(:,i); % initial condition for this time step
    u = u_test(:,i); % Assume constant u at average of time interval
%     [t_1,x_1] = ode45(@(t_1,x_1) cartpend(x_1,u), t_test(i:i+1), x00);
    [t_1,x_1] = ode45(@(t_1,x_1) (A_lin*x_1 + B_lin*(u - u_bar)), t_test(i:i+1), x00);
    x_hat3(:,i+1) = x_1(end,:);
end

y_hat3 = x_hat3([1,3],:); % Extract measurements

MAE_lin = sum(abs(y_hat3 - y_test), 2)./N_test % MAE for linearised model

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

plot(t, u_data, ':', 'LineWidth', 1);
plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_train, y_hat2, '--', 'LineWidth', 1); % Plot only non-delay coordinate  
plot(t_test, y_hat3, ':', 'LineWidth', 1); % Plot linearised model results
plot((D + t(N-N_test-N_train)).*[1,1], ylim, 'r');
plot(t(N-N_test-N_train).*[1,1], ylim, 'k');
plot(t(N-N_test).*[1,1], ylim, 'k');
title('Training and Testing data vs Model');
% legend('x', 'theta', 'input', 'x_hat', 'theta_hat', 'D', 't(final sample)')
hold off;

toc(total_timer);

disp('-------------------')
disp('END of HAVOK script')


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

function dx = cartpend(x,u)
    
    %CARTPEND Models a continuous system of a pendulem on a cart.
    %   based on Steve Brunton code. See youtube.com/watch?v=qjhAAQexzLg&list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m&index=12
    %   x  = state vector [x; x_dot; theta; theta_dot]
    %   dx = derivative of state vector
    %   u  = input vector [f]
    %   m  = mass of pendulem end
    %   M  = mass of cart
    %   L  = length of pendulem rod
    %   g  = acceleration due to gravity
    %   d  = damping coef of friction on cart
    %,m,M,L,g,d,u
    
    u_disturb = 5;
    
    m = 2;
    M = 4;
    L = 1;
    g = -9.81;
    d = 5;
    
    % Derivatives
    dx = zeros(4,1);

    Sx = sin(x(3));
    Cx = cos(x(3));
    D = m*L*L*(M+m*(1-Cx^2));

    % Equations from derive_cartpend.m
    dx(1,1) = x(2);
    dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*(u + u_disturb);
    dx(3,1) = x(4);
    dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*(u + u_disturb); % +.01*randn;

end




