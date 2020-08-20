%% HAVOK with control - of floating pendulum 2D
% Estimate linear model from data
% Partial state feedback

<<<<<<< HEAD
%% Try to save previous results of random search and continue with them

%% Read data

close all;

total_timer = tic;

% load('floating_pend_2D_random_1.mat') % Load simulation data
% load('floating_pend_2D_PI_z_control_and_random_1.mat') % simulation data with z controlled by PI controller

u_data  = out.u.Data';
x_data  = out.x.Data';
% NB! set which variables are measured
y_data  = x_data([1:3],:); % Measurement data (x,z, theta)
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
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

%% Parameters
% Very dependant on choice of p, r, q

sigma = 0.001; % Noise standard deviation
N_train = 3000; % Number of sampels in training data
c = 1; % Column spacing of Hankel matrix (for multiscale dynamics)
d = 1; % Row spacing of Hankel matrix (for multiscale dynamics)
% w; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
% p % Truncated rank of system of Omega
% q % number of delays
% r = Truncated rank of system of X2

% Add noise
rng('default'); % Initialise random number generator
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));

%% Load saved results
estimation_name = 'HAVOK'; % Name of prediction model
plant_name = 'floating_pend_2D'; % Name of prediction model
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
save_file = ['Data/', plant_name, '_', estimation_name, '_N_train_vs_error', '_sig=', sig_str, '.mat'];

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
    error(['Saved results file:', newline, save_file, newline, 'does not exist'])  
end

% r = p - 2;
r = p - l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)
w = N_train - q; % num columns of Hankel matrix
D = (q-1)*d*Ts; % Delay duration (Dynamics in delay embedding)

% Training data - Last sample of training is first sample of testing
y_train = y_data_noise(:, end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:, end-N_test-N_train+2:end-N_test+1);
t_train = t(:, end-N_test-N_train+2:end-N_test+1);

timer_model = tic; % Start timer for this model evaluation

% Step 1: Collect and construct the snapshot matrices:
% According to: Discovery of Nonlinear Multiscale Systems: Sampling Strategies and Embeddings
% pg. 15 Delay spacing for multiscale dynamics
=======
%% Read data

load('floating_pend_2D_random_1.mat') % Load simulation data
u_data  = out.u.Data';
x_data  = out.x.Data';

n = size(x_data,1); % number of states
I = eye(n);
C = I([1,2,3],:); % Measurement matrix (Only measures x,z,theta)
y_data  = C*x_data; % Measure x and theta
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
t  = out.tout';
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% [X_p,Y_delays] = meshgrid(1:delays(end), 1:delays(end)); % Values for surface plot
% RMSE_matrix = zeros(delays(end), delays(end)); % Empty matrix of errors

% Very dependant on choice of delays, p, r, q
p = 140; % Truncated rank of system
c = 1; % Column spacing of Hankel matrix
d = 1; % Row spacing of Hankel matrix
q = 1500; % number of delays
w = 4000; % (named 'p' in Multiscale paper) number of columns in Hankel matrix
sigma = 0; % Noise standard deviation

numel_H = q*m*w;
log10(numel_H)
time_predict = 6e-6*numel_H

final_sample = (q-1)*d + (w-1)*c + 2 % Last sample used in Hankel matrix
assert(final_sample <= N, 'Not enough data. Change q, d, w, or c'); % otherwise index out of bounds 
D = (q-1)*d*Ts; % Delay duration

% Train/Test split
N_train = final_sample; % Num of data samples for training, rest for testing
y_train = y_data(:,1:N_train);
u_train = u_data(:,1:N_train);
t_train = t(:,1:N_train);

N_test = N - N_train + 1; % Num of data samples for testing
y_test = y_data(:,N_train:end); % One sample overlaps for initial condition
u_test = u_data(:,N_train:end);
t_test = t(:,N_train:end);

% Add noise
y_train = y_train + sigma*randn(size(y_train));

%%

for q = q % Loop through delays
    for p = p % Loop truncateed rank

r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)
tic;

% Step 1: Collect and construct the snapshot matrices:
% According to: Discovery of Nonlinear Multiscale Systems: Sampling Strategies and Embeddings
% pg. 15
% Delay spacing for multiscale dynamics
disp(1)
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd
X = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
X2 = zeros(q*m,w); % X one step into future
for row = 0:q-1 % Add delay coordinates
    X(row*m+1:(row+1)*m, :) = y_train(:, row*d + (0:w-1)*c + 1);
    X2(row*m+1:(row+1)*m, :) = y_train(:, row*d + (0:w-1)*c + 2);
end

Upsilon = u_train(:, row*d + (0:w-1)*c + 1); % Upsilon, same indexes as last X row

Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

<<<<<<< HEAD
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
=======
% Step 2: Compute and truncate the SVD of the input space Omega
disp(2)
[U,S,V] = svd(Omega, 'econ');
figure, semilogy(diag(S), 'x'), hold on;
title('Singular values of Omega, showing p truncation')
plot(p,S(p,p), 'ro'), hold off;
figure
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd
% for i=1:20    
%     plot(V(:,i));
%     pause
% end
<<<<<<< HEAD

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
=======
U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
S_tilde = S(1:p, 1:p);
V_tilde = V(:, 1:p);
U1_tilde = U_tilde(1:q*m, :);
U2_tilde = U_tilde(q*m+1:end, :);

% Step 3: Compute the SVD of the output space X'
disp(3)
[U,S,V] = svd(X2, 'econ');
figure, semilogy(diag(S), 'x'), hold on;
title('Singular values of X2, showing r truncation')
plot(r,S(r,r), 'ro'), hold off;
U_hat = U(:, 1:r); % Truncate SVD matrixes of X2
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);

% Step 4: Compute the approximation of the operators G = [A B]
disp(4)
A_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U1_tilde'*U_hat;
B_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U2_tilde';
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd

% x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
% x = U_hat*x_tilde, Transform to original coordinates
% x_tilde = U_hat'*x, Transform to reduced order coordinates
% Here x is augmented state

A = U_hat*A_tilde*U_hat';
B = U_hat*B_tilde;

<<<<<<< HEAD
if (sum(abs(eig(A)) > 1) ~= 0) % If eigenvalues are unstable
    error('eigenvalues are unstable');
end

% Time taken to train this model
time = toc(timer_model) % Add time taken in q loop before p chosen

% x_augmented(k+1) = A*x_aug(k) + B*u(k)

%% Compare to testing data
=======
if (sum(abs(eig(A)) > 1) ~= 0) % If some eigenvalues are unstable due to machine tolerance
    disp('Unstable eigenvalues')
end

% If some eigenvalues are unstable due to machine tolerance,
% Scale them to be stable
count = 0;
while (sum(abs(eig(A)) > 1) ~= 0) 
    count = count+1
    [Ve,De] = eig(A);
    unstable = abs(De)>1; % indexes of unstable eigenvalues
    De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-16+count); % Normalize all unstable eigenvalues (set abs(eig) = 1)
    A_tilde = Ve*De*inv(Ve); % New A with margininally stable eigenvalues
    A_old = A;
    A = real(A);
    if(count>10)
        'break'
        break
    end
end
assert(sum(abs(eig(A)) > 1) == 0, 'Unstable eigenvalues'); % Check if all eigenvalues are stable (magnitude <= 1)

disp('Model computed')
toc;
% x_augmented(k+1) = A*x_aug(k) + B*u(k)

% Ignore eigenmodes Step 5 and 6

%% Compare to testing data
disp(5)
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd
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

<<<<<<< HEAD
% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state


%% Compare to training data
disp('Compare to training data')
% Initial conditions
y_hat_02 = zeros(q*m,1);
for row = 0:q-1 % Create first column of spaced Hankel matrix
    y_hat_02(row*m+1:(row+1)*m, 1) = y_train(:, row*d + 1);
=======
toc;

%% Compare to training data
disp(6)
% Initial conditions
y_hat_02 = zeros(q*m,1);
for row = 0:q-1 % Create first column of spaced Hankel matrix
    y_hat_02(row*m+1:(row+1)*m, 1) = y_data(:, row*d + 1);
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd
end
k_start = row*d + 1; % First k to start at

Y_hat2 = zeros(length(y_hat_0),N_train); % ??? Estimated X from model
Y_hat2(:,k_start) = y_hat_02; % Initial conditions, insert at first k
for k = k_start:N_train-1
    Y_hat2(:,k+1) = A*Y_hat2(:,k) + B*u_train(:,k);
end
y_hat2 = Y_hat2(end-m+1:end, :); % Extract only non-delay time series (last m rows)

<<<<<<< HEAD
disp('Run model on training data')

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

% plot(t, u_data, ':', 'LineWidth', 1);
plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_train, y_hat2, '--', 'LineWidth', 1); % Plot only non-delay coordinate  
plot((D + t(N-N_test-N_train)).*[1,1], ylim, 'r');
plot(t(N-N_test-N_train).*[1,1], ylim, 'k');
plot(t(N-N_test).*[1,1], ylim, 'k');
=======
%% Plot data vs model
figure;
plot(t, y_data); 
hold on;
plot(t, u_data, ':', 'LineWidth', 1);
plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_train, y_hat2, '--', 'LineWidth', 1); % Plot only non-delay coordinate  
plot([D D], ylim, 'r');
plot([t(N_train) t(N_train)], ylim, 'k');
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd
title('Training and Testing data vs Model');
% legend('x', 'theta', 'input', 'x_hat', 'theta_hat', 'D', 't(final sample)')
hold off;

toc(total_timer);

disp('-------------------')
disp('END of HAVOK script')

toc;
    end
end

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

<<<<<<< HEAD
function new_array = insert(array, index, entry)
    if index == 1 % to avoid index-1 = 0
        new_array = [entry, array];
    else
        new_array = [array(:, 1:index-1), entry, array(:, index:end)];
    end
end


=======
>>>>>>> 635cc56f06e099eba84aca5266d04b9d329375cd




