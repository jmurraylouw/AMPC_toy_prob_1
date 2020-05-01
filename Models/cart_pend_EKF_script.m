% Applies Extended Kalman Filter for 
% simueltaneous state and parameter estimation
% on data from cart_pend

% Read simulation data
x_data = out.x.Data';
y_data = out.y.Data';
u_data = out.u.Data';
t = out.x.Time';
Ts = t(2)-t(1);

% Dimensions
[nx, n_time] = size(x_data)
[ny, n_time] = size(y_data)
[nu, n_time] = size(u_data)

% Initialise
x0 = [0; 0; 0; 0.5; 2; 10; 10; 4];
nx = length(x0);
P0 = 0.5*eye(nx);
u0 = 0;
x_hat = x0;
P = P0;
u = u0;

Q = 0.00001*eye(nx); % Model uncertainty
R = 0.0001*eye(ny); % Measurement uncertainty

% System definition
m = 1;
M = 5;
L = 2;
g = -9.81;
d = 1;

f = @cartpend; % Function handle
g = @measure; % Measurement function handle
A = jaccsd(f,x_hat,u);
F = A;
B = [0; 1/M; 0; -1/(M*L); 0]; % ?? Later change jaccsd to also compute B
C = jaccsd(g,x_hat,u);
H = C;
D = 0;

% Extrapolate
x_hat_dwork = x_hat + f(x_hat,u)*Ts; % Numeric integration to extrapolate state

F = jaccsd(f,x_hat,u); % Calculate Jacobian of continuous system
Phi = eye(nx) + F*Ts + 1/2*(F*Ts)^2; % ??? where is this from? 2nd order Taylor expansion? (continuous to discrete)
P_dwork = Phi*P*Phi' + Q; % Extrapolate uncertainty

x_hat_data = zeros(nx, n_time); % Assign memory beforehand

% Apply EKF at every timestep
for n = 1:1:n_time-1
    % Measurement
    y = y_data(n);
    u = u_data(n);
    
    % Get saved data
    x_hat = x_hat_dwork;
    P = P_dwork;
    
    % Update
    H = jaccsd(g,x_hat,0); % Linearise measurement function
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

plot_rows = [1 3];
figure
plot(t, x_data(plot_rows,:)); 
hold on
plot(t, x_hat_data(plot_rows,:));
% plot(t, y_data);
hold off;
legend('Actual x', 'Actual theta', 'Estimate x','Estimate theta')


function dx = cartpend(x,u)
% Adapted from code by Steve Brunton
% x contains state, input and parameters
% x = [x;
%     x_dot;
%     theta;
%     theta_dot;
%     m;]

% Parameters
m = x(6); % 1
M = x(8); % 5
L = x(5); % 2
g = -9.81;
d = x(7); % 10

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
y(1) = x(1);
y(2) = x(3);
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
% ?? Maybe later add calculation of B also
h = n*eps;
for k=1:n
    x1 = x;
    x1(k) = x1(k)+ h*1i;
    J(:,k) = imag(f(x1,u))/h;
end
end



