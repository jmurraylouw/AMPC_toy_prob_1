%% HAVOK with control - of cart pendulum
% Estimate System matrixes with moving window of data in real time
% Partial state feedback

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
plot(t,x_data)
%%
% Validation data
load('cartpend_data_4.mat') % Load validation data
u_valid  = out.u.Data';
x_valid  = out.x.Data';
I = eye(n);
y_valid  = I([1, 3],:)*x_valid; % Measure x and theta
t_valid  = out.tout';
Ts_valid = t_valid(2)-t_valid(1); % Sample time of data

% Very dependant on choice of delays, p, r

samples = floor(N*0.6);
% delays = 1; % Number of delay cordinates, including y_data(1:samples+1)
% p = 1; % Reduced rank of Omega svd
% [X_p,Y_delays] = meshgrid(1:delays(end), 1:delays(end)); % Values for surface plot
% RMSE_matrix = zeros(delays(end), delays(end)); % Empty matrix of errors

c = 1; % Column spacing of Hankel matrix
d = 1; % Row spacing of Hankel matrix

q = 4; % number of delays
w = samples; % number of columns in Hankel matrix

for q = 1000 % Loop through delays
    for p = 7

r = p-l; % Reduced rank of X2 svd, r < p, (minus number of inputs from rank)
tic;
% Noticed error increased for increasing tau
% i.e Y = y(k)      y(k+1)      y(k+2)...
%         y(k+tau)  y(k+tau+1)  y(k+tau+2)...

assert(samples+q <= N); % otherwise index out of bounds 

% Step 1: Collect and construct the snapshot matrices:
1
% According to: Discovery of Nonlinear Multiscale Systems: Sampling Strategies and Embeddings
% pg. 15
% Delay spacing for multiscale dynamics

X = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
X2 = zeros(q*m,w); % X one step into future
for row = 0:q-1 % Add delay coordinates
    X(row*m+1:(row+1)*m, :) = y_data(:, row*d + (0:w-1)*c + 1);
    X2(row*m+1:(row+1)*m, :) = y_data(:, row*d + (0:w-1)*c + 2);
end

Upsilon = u_data(:, row*d + (0:w-1)*c + 1); % Upsilon, same indexes as last X row

Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

% Step 2: Compute and truncate the SVD of the input space Omega
2
[U,S,V] = svd(Omega, 'econ');
figure(1), semilogy(diag(S), 'x')
U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
S_tilde = S(1:p, 1:p);
V_tilde = V(:, 1:p);
U1_tilde = U_tilde(1:q*m, :);
U2_tilde = U_tilde(q*m+1:end, :);

% Step 3: Compute the SVD of the output space X'
3
[U,S,V] = svd(X2, 'econ');
figure(2), semilogy(diag(S), 'x')
U_hat = U(:, 1:r); % Truncate SVD matrixes of X2
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);

% Step 4: Compute the approximation of the operators G = [A B]
4
A_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U1_tilde'*U_hat;
B_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U2_tilde';

% x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
% x = U_hat*x_tilde, Transform to original coordinates
% x_tilde = U_hat'*x, Transform to reduced order coordinates
% Here x is augmented state

A = U_hat*A_tilde*U_hat';
B = U_hat*B_tilde;

% count = 0;
% while (sum(abs(eig(A)) > 1) ~= 0) % If some eigenvalues are unstable due to machine tolerance
%     count=count+1;
%     [Ve,De] = eig(A);
%     unstable = abs(De)>1; % indexes of unstable eigenvalues
%     De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-10+count); % Normalize all unstable eigenvalues (set abs(eig) = 1)
%     A = Ve*De*inv(Ve); % New A with margininally stable eigenvalues
%     if(count>10)
%         'break'
%         break
%     end
% end
% assert(sum(abs(eig(A)) > 1) == 0, 'Unstable eigenvalues'); % Check if all eigenvalues are stable (magnitude <= 1)

% x_augmented(k+1) = A*x_aug(k) + B*u(k)

% Ignore eigenmodes Step 5 and 6

%% Compare to validation data
5
% Initial conditions
x_hat_0 = zeros(q*m,1);
for row = 0:q-1
    x_hat_0(row*m+1:(row+1)*m, 1) = y_valid(:, row*d + 1); % Add first row of spaced Hankel matrix
end
k_start = row*d + 1; % First k to start at

X_hat = zeros(length(x_hat_0),N); % ??? Estimated X from model
X_hat(:,k) = x_hat_0; % Initial conditions
for k = (row*d + 1):N-1
    X_hat(:,k+1) = A*X_hat(:,k) + B*U_data(k);
end

toc;
figure(3), plot(t_valid, x_hat([1, 2],:), '--'); % Plot only non-delay coordinate
hold on;
plot(t_valid, x_valid([1, 3],:))

hold off;

% figure(4), plot(t_cut, x_data_cut([1, 3],:)-x_hat([1, 2],:), 'r');
% title('Error')
p
RMSE1 = mean(abs(x_valid(1,:).^2 - x_hat(1,:).^2))
RMSE_matrix(q,r) = RMSE1;


    end
end

% surf(X_p,Y_delays,abs(RMSE_matrix))
% min_error = min(min(RMSE_matrix(RMSE_matrix~=0)))
% 
% [col_mins, row_indexes] = min(RMSE_matrix(RMSE_matrix~=0));
% [abs_min, col_index] = min(col_mins);
% delays = row_indexes(col_index)
% p = col_index
% %%
% RMSE_matrix(RMSE_matrix>2) = RMSE_matrix(RMSE_matrix>2)*0;
% bar3(RMSE_matrix)
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





