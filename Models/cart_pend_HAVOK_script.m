%% HAVOK with control - of cart pendulum
% Estimate System matrixes with moving window of data in real time
% Partial state feedback

%% Variables for simulation
x0 = [1; 0; -0.2; -0.8];

%% Read data
% Training data
load('cartpend_data_3.mat') % Load training data
u_data  = out.u.Data';
x_data  = out.x.Data';
n = size(x_data)*[1; 0]; % number of states
I = eye(n);
y_data  = I([1, 3],:)*x_data; % Measure x and theta
n = size(x_data)*[1; 0]; % number of states
m = size(y_data)*[1; 0]; % number of measurements
l = size(u_data)*[1; 0]; % number of inputs
t  = out.tout';
Ts = t(2)-t(1);     % Sample time of data
N  = length(t); % Number of data samples
% plot(t,x_data)

% Validation data
load('cartpend_data_4.mat') % Load validation data
u_valid  = out.u.Data';
x_valid  = out.x.Data';
I = eye(n);
y_valid  = I([1, 3],:)*x_valid; % Measure x and theta
t_valid  = out.tout';
Ts_valid = t_valid(2)-t_valid(1);     % Sample time of data

% Very dependant on choice of delays, p, r

samples = floor(N*0.8);
delays = 20; % Number of delay cordinates, including y_data(1:samples+1)
p = 9; % Reduced rank of Omega svd
r = p-1; % Reduced rank of X2 svd, r < p
tic;
tau = 1; % Sample number shift of delay cordinates.
% Noticed error increased for increasing tau
% i.e Y = y(k)      y(k+1)      y(k+2)...
%         y(k+tau)  y(k+tau+1)  y(k+tau+2)...

assert(samples+delays*tau <= N); %% otherwise index out of bounds 

% Step 1: Collect and construct the snapshot matrices:
X = []; % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
for i = 1:tau:delays*tau
    X = [y_data(:, i:samples+i); X];  
end

na = size(X)*[1; 0]; % Length of augmented state vector

X2 = X(:, 2:end); % X advanced 1 step into future
X = X(:, 1:end-1); % Cut off last sample

Upsilon = u_data(:, (delays*tau):(samples + delays*tau - 1)); % Upsilon

Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

% Step 2: Compute and truncate the SVD of the input space Omega
[U,S,V] = svd(Omega, 'econ');
% figure(1), plot(diag(S))

% plot(diag(S), 'o')
U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
S_tilde = S(1:p, 1:p);
V_tilde = V(:, 1:p);
U1_tilde = U_tilde(1:na, :);
U2_tilde = U_tilde(na+1:end, :);

% Step 3: Compute the SVD of the output space X'
[U,S,V] = svd(X2, 'econ');
% figure(2), plot(diag(S))

% plot(diag(S), 'o')
U_hat = U(:, 1:r); % Truncate SVD matrixes of X2
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);

% Step 4: Compute the approximation of the operators G = [A B]
A_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U1_tilde'*U_hat;
B_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U2_tilde';

% x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
% x = U_hat*x_tilde, Transform to original coordinates
% x_tilde = U_hat'*x, Transform to reduced order coordinates
% Here x is augmented state

A = U_hat*A_tilde*U_hat';
B = U_hat*B_tilde;

assert(sum(abs(eig(A)) > 1) == 0); % Check if all eigenvalues are stable (magnitude <= 1)

% x_augmented(k+1) = A*x_aug(k) + B*u(k)

% Ignore eigenmodes Step 5 and 6

%% Compare to validation data

% Initial conditions
x_hat_0 = [];
for i = 1:tau:delays*tau
    x_hat_0 = [y_valid(:,i); x_hat_0];
end

t_valid_cut = t_valid(delays*tau:end); % Cut off first part to match model
x_valid_cut = x_data(:,delays*tau:end);


% x_tilde_0 = U_hat'*x_hat_0; % Convert to tilde coordinate system
% x_tilde_hat = run_model(A_tilde,B_tilde,u_data,t_cut,x_tilde_0);
% x_hat = U_hat*x_tilde_hat; % Convert back to original coordinate system

x_hat = run_model(A,B,u_data,t_valid_cut,x_hat_0); 
toc;
figure(3), plot(t_valid_cut, x_hat([1, 2],:), '--'); % Plot only non-delay coordinate
hold on;
plot(t_valid_cut, x_valid_cut([1, 3],:))
hold off;

% figure(4), plot(t_cut, x_data_cut([1, 3],:)-x_hat([1, 2],:), 'r');
% title('Error')

delays
tau
MSE_dmdc = mean(((x_data_cut(1,:) - x_hat(1,:)).^2)')'

%%
stop


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
    X_hat(:,1) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
end





