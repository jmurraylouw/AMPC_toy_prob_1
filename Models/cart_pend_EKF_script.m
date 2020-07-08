% Applies Extended Kalman Filter for 
% simueltaneous state and parameter estimation
% on data from cart_pend

% System definition
f = @cartpend; % Function handle
g = @measure; % Measurement function handle

% Read simulation data
load('cartpend_random_1');
x_data = out.x.Data';
y_data = g(x_data);
u_data = out.u.Data';
t = out.x.Time';
Ts = t(2)-t(1);
N = length(t);

% Train/Test split
N_train = 3000; % Num of data samples for training, rest for testing
y_train = y_data(:,1:N_train);
u_train = u_data(:,1:N_train);
t_train = t(:,1:N_train);

N_test = N - N_train + 1; % Num of data samples for testing
x_test = x_data(:,N_train:end); % One sample of testing data overlaps for initial condition
y_test = y_data(:,N_train:end); 
u_test = u_data(:,N_train:end);
t_test = t(:,N_train:end);

% Add noise
sigma = 0.0;
y_train = y_train + sigma*randn(size(y_train));

% Initialise
% x = [x, x_dot, theta, theta_dot, L, m, d]
x0 = [0; 0; 0; 0; 2; 4; 1; 5]; % Correct initial values
x0 = x0 + 0*randn(size(x0)); % Add random alue to initial guesses
nx = length(x0); % 4 states, 3 paramters
ny = length(g(x0)); % x and theta
u0 = 0;
nu = length(u0);

P0 = 0.1*eye(nx); % Initial guess uncertainty
Q = diag([0; 0.00001; 0; 0.00001; 0.00001; 0.00001; 0.00001; 0.00001]); % Model uncertainty
R = sigma*eye(ny); % Measurement uncertainty

x_hat = x0;
P = P0;
u = u0;

% Extrapolate
x_hat_dwork = x_hat + f(x_hat,u)*Ts; % Numeric integration to extrapolate state

F = jaccsd(f,x_hat,u); % Calculate Jacobian of continuous system
Phi = eye(nx) + F*Ts + 1/2*(F*Ts)^2; % ??? where is this from? 2nd order Taylor expansion? (continuous to discrete)
P_dwork = Phi*P*Phi' + Q; % Extrapolate uncertainty

x_hat_data = zeros(nx, N_train); % Assign memory beforehand

% Apply EKF at every timestep
for n = 1:1:N_train
    % Measurement
    y = y_train(:,n);
    u = u_train(:,n);
    
    % Get saved data
    x_hat = x_hat_dwork;
    P = P_dwork;
    
    % Update
    H = jaccsd(g,x_hat,u); % Linearise measurement function
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Update estimate with measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,n) = x_hat;
    
    % Extrapolate for next time step
    x_hat = x_hat + f(x_hat,u)*Ts; % Numeric integration (extrapolate state)

    F = jaccsd(f,x_hat,0); % Calculate Jacobian of continuous system
    Phi = eye(nx) + F*Ts + 0.5*(F*Ts)^2; % 2nd order Taylor expansion (continuous to discrete)
    P = Phi*P*Phi' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end

% Actual parameters
m = 2; % Mass of payload
M = 4; % Mass of cart
L = 1; % Length of pendulum
d = 5; % Linear damping coef on cart

% Plot parameter estimation
param_rows = [5:8];
figure
plot(t_train, m*ones(1,N_train)); hold on
plot(t_train, M*ones(1,N_train));
plot(t_train, L*ones(1,N_train));
plot(t_train, d*ones(1,N_train));
plot(t_train, x_hat_data(param_rows,:));
hold off;
legend('m','M','L','d','m est','M est','L est','d est');
title('Parameter estimation');

% Extract final model parameters
m = x_hat_data(5,end);
M = x_hat_data(6,end);
L = x_hat_data(7,end);
d = x_hat_data(8,end);

%% Testing
% Run model on unseen testing data and compare to actual measurements

x0 = [x_test(:,1);m;M;L;d]; % Initial condition
x_hat = zeros(nx,N_test); % Empty prediction matrix
x_hat(:,1) = x0; 
% Generate data with SINDY-PI model
% Solve for small intervals with constant u
for i=1:N_test-1
    x0 = x_hat(:,i); % initial condition for this time step
    u = u_test(:,i); %(U_test(i,:) + U_test(i+1,:))/2; % Assume constant u at average of time interval
    [t_1,x_1] = ode45(@(t_1,x_1) cartpend(x_1,u), t_test(i:i+1), x0);
    x_hat(:,i+1) = x_1(end,:);
end

y_hat = x_hat([1,3],:);
y_test = x_test([1,3],:);

% Vector of Root Mean Squared Error on testing data
RMSE = sqrt(sum((y_hat - y_test).^2, 2)./N_test) % Each row represents RMSE for measured state [RMSE_x]

state_rows = [1,3];
figure
plot(t, y_data); hold on
plot(t_train, x_hat_data(state_rows,:), '--'); % Plot state estimation on trainging data
plot(t_test, y_hat, '--');
plot([t(N_train) t(N_train)], ylim, 'k');
hold off;
title('State estimation');



function dx = cartpend(x,u)
% Adapted from code by Steve Brunton
% x contains state, input and parameters
% x = [x;
%     x_dot;
%     theta;
%     theta_dot;
%     L;]

% Parameters
m = x(5);% 2;
M = x(6); % 4;
L = x(7); % 1
g = -9.81;
d = x(8); % 5;

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

nx = length(x);

dx = zeros(nx,1); % Assign memory space
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;
end

function y = measure(x,u)
    % Measurement function    
    y(1,:) = x(1,:);
    y(2,:) = x(3,:);
end


function J=jaccsd(f,x,u) % ??? Maybe should use simbolic diff for more exact
% JACCSD Jacobian through complex step differentiation
% By Yi Cao at Cranfield University, 02/01/2008
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
f_x = f(x,u);
n = numel(x);
m = numel(f_x);
J = zeros(m,n);
h = n*eps;
for k=1:n
    x1 = x;
    x1(k) = x1(k)+ h*1i;
    J(:,k) = imag(f(x1,u))/h;
end
end



