% KF
%{
x_hat_data = zeros(nx, n_time);

sigma_a = 0.1; % Std dev of acceleration/force applied to model
Q = G*(sigma_a^2)*G'
R = 0.01;

% Initialise
x0 = [0; 0];
P0 = [0, 0; 0, 0];
x_hat = x0;
P = P0;

x_hat_dwork = F*x_hat + G*u; % Extrapolate state
P_dwork = F*P*F' + Q; % Extrapolate uncertainty

 
for n = 1:1:n_time-1
    % Measurement
    y = y_data(n);
    u = u_data(:,n);
    
    % Get saved data
    x_hat = x_hat_dwork;
    P = P_dwork;
    
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Update estimate with measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,n) = x_hat;
    
    x_hat = F*x_hat + G*u; % Extrapolate state
    P = F*P*F' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end

plot(t, x_data(1,:)); hold on;

plot(t, x_hat_data(1,:));
%plot(t, y_data);
hold off;
legend('Actual', 'Estimate', 'Measured')
%}
f = @cartpend; % Function handle
x = [0; 0; 0; 0];

m = 1;
M = 5;
L = 2;
g = -10;
d = 1;
c = [m; M; L; g; d];

u = 1;
[f_x,J] = jaccsd(f,x,c,u)

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

function dx = cartpend(x,c,u)
% Adapted from code by Steve Brunton
% Parameter vector, c
m = c(1);
M = c(2);
L = c(3);
g = c(4);
d = c(5);

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

dx = zeros(4,1); % Assign memory space
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;
end

function [f_x,J]=jaccsd(f,x,c,u)
% JACCSD Jacobian through complex step differentiation
% By Yi Cao at Cranfield University, 02/01/2008
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
f_x = f(x,c,u);
n = numel(x);
m = numel(f_x);
J = zeros(m,n);
h = n*eps;
for k=1:n
    x1 = x;
    x1(k) = x1(k )+ h*1i;
    J(:,k) = imag(f(x1,c,u))/h;
end
end



