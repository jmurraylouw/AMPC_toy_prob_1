close all;

% Read simulation data
x_data = out.x.Data';
y_data = out.y.Data';
u_data = out.u.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data)
[ny, n_time] = size(y_data)
[nu, n_time] = size(u_data)

% Initialise
x0 = [0; 0; 0; 1];
P0 = 0.01*eye(nx);
x_hat = x0;
P = P0;

% System definition
m = 1;
M = 5;
L = 2;
g = -9.81;
d = 1;

f = @cartpend; % Function handle
[f_x,A] = jaccsd(f,x0,u);
B = [0; 1/M; 0; s*1/(M*L); 0]; % ?? Later change jaccsd to also compute B
C = eye(nx);
% C = [1 0 0 0 0 0;
%      0 0 1 0 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1]; % Measure x, theta, u, m
D = 0;

A

sys_c = ss(A,B,C,D); % Continuous system
sys_d = c2d(sys_c, Ts, 'zoh'); % Discrete system
[F,B,H,D] = ssdata(sys_d);

sigma_a = 0.1; % Std dev of acceleration/force applied to model
Q = 0.01*eye(nx); % Model uncertainty
R = 0.0001*eye(ny); % Measurement uncertainty

% Extrapolate
x_hat_dwork = F*x_hat + B*u; % Extrapolate state
P_dwork = F*P*F' + Q; % Extrapolate uncertainty

x_hat_data = zeros(nx, n_time); % Assign memory beforehand

% Apply EKF at every timestep
for n = 1:1:n_time-1
    % Measurement
    y = y_data(n);
    
    % Get saved data
    x_hat = x_hat_dwork;
    x_hat_dwork;
    P = P_dwork;
    
    % Update
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Update estimate with measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,n) = x_hat;
    
    % Linearise system equations
    [f_x,A] = jaccsd(f,x_hat,u);
    sys_c = ss(A,B,C,D); % Continuous system
    sys_d = c2d(sys_c, Ts, 'zoh'); % Discrete system
    [F,B,H,D] = ssdata(sys_d);
    
    % Extrapolate for next time step
    x_hat = F*x_hat + B*u; % Extrapolate state
    P = F*P*F' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end

plot_rows = [1 3 5];
figure
plot(t, x_data(plot_rows,:)); hold on

% plot(t, y_data);

plot(t, x_hat_data(plot_rows,:));

hold off;
legend('Actual', 'Estimate', 'Measured')


function dx = nl_msd(x,c,u)
% NON-LINEAR MASS SPRING DAMPER
% INPUT:
% x = state vector [x, x_dot]
% u = input vector [f]
% c = paramter vecotr [m, b, k]
% 
% OUTPUT:
% dx = derivative of state vector [x_dot, x_dotdot]

m = c(1); % Mass
b = c(2); % Damper coef
k = c(3); % Coef of non-linear spring

dx = zeros(2,1); % Assign memory
dx(1,1) = x(2);
dx(2,1) = 1/m*(-k*x(1)^3 - b*x(2) + u);
end


function dx = cartpend(x,u)
% Adapted from code by Steve Brunton
% x contains state, input and parameters
% x = [x;
%     x_dot;
%     theta;
%     theta_dot;
%     m;]

% Parameter vector
m = x(5);
M = 5;
L = 2;
g = -9.81;
d = 10;

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

dx = zeros(4,1); % Assign memory space
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;
dx(5,1) = 0;

end
% Old cartpend:
%{
function dx = cartpend(x)
% Adapted from code by Steve Brunton
% Parameter vector, c
m = c(1);
M = c(2);
L = c(3);
g = -9.81;
d = c(4);

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

dx = zeros(4,1); % Assign memory space
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;
end
%}

function [f_x,J]=jaccsd(f,x,u) % ??? Maybe should use simbolic diff for more exact
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



