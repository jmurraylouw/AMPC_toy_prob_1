% Sets up workspace for cart_pend_EKF SIMULINK model to run.
% Applies Extended Kalman Filter for 
% simueltaneous state and parameter estimation

Ts = 0.01;

% Dimensions
nx = 4 + 3; % 4 states, 3 paramters
ny = 2; % x and theta
nu = 1;

% Initialise
% x = [x, x_dot, theta, theta_dot, L, m, d, M]
x0 = [0; 0; 0; 0.5; 2; 10; 10; 4];
nx = length(x0);
P0 = 0.5*eye(nx);
u0 = 0;

Q = 0.00001*eye(nx); % Model uncertainty
R = 0.0001*eye(ny); % Measurement uncertainty

% Function handles
f = @cartpend; % Function handle
g = @measure; % Measurement function handle

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



