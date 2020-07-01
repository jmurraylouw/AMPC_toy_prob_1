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
y_test = y_data(:,N_train:end); % One sample overlaps for initial condition
u_test = u_data(:,N_train:end);
t_test = t(:,N_train:end);

% Add noise
sigma = 0.1;
y_train = y_train + sigma*randn(size(y_train));

% Actual parameters
m = 2; % Mass of payload
M = 4; % Mass of cart
L = 1; % Length of pendulum
d = 5; % Linear damping coef on cart

% Dimensions
nx = size(x_data,1); % number of states
ny = size(y_data,1); % number of measurements
nu = size(u_data,1); % number of inputs

% Initialise
x0 = [1; -0.2; -0.5; 0.8; 1.1];
nx = length(x0);
P0 = 0.5*eye(nx);
u0 = 0;
x_hat = x0;
P = P0;
u = u0;

% Uncertainty values
Q = diag([0; 0; 0.000001; 0.0000001; 0.001]); % Model uncertainty
R = 0.1*eye(ny); % Measurement uncertainty

% Extrapolate
x_hat_dwork = x_hat + f(x_hat,u)*Ts; % Numeric integration to extrapolate state

F = jaccsd(f,x_hat,u); % Calculate Jacobian of continuous system
Phi = eye(nx) + F*Ts + 1/2*(F*Ts)^2; % ??? where is this from? 2nd order Taylor expansion? (continuous to discrete)
P_dwork = Phi*P*Phi' + Q; % Extrapolate uncertainty

x_hat_data = zeros(nx, N_train); % Assign memory beforehand

% Apply EKF at every timestep
for n = 1:1:N_train-1
    % Measurement
    y = y_train(n);
    u = u_train(n);
    
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

plot_rows = [1,2,3,4];
figure
plot(t, y_data); hold on
plot(t, m*ones(1,N));
plot(t_train, x_hat_data(plot_rows,:));
hold off;
legend('Actual x', 'Actual theta', 'Estimate x','Estimate theta')


function dx = cartpend(x,u)
% Adapted from code by Steve Brunton
% x contains state, input and parameters
% x = [x;
%     x_dot;
%     theta;
%     theta_dot;
%     L;]

% Parameters
m = 2; %x(6);
M = 4; %x(8);
L = x(5);
g = -9.81;
d = 5; %x(7);

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
    y(2,:) = x(2,:);
    y(3,:) = x(3,:);
    y(4,:) = x(4,:);
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



