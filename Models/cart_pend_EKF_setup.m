% Sets up workspace for cart_pend_EKF SIMULINK model to run.
% Applies Extended Kalman Filter for 
% simueltaneous state and parameter estimation

Ts = 0.01;

% Function handles
f = @cartpend; % Function handle
g = @measure; % Measurement function handle

% Initialise
% x = [x, x_dot, theta, theta_dot, L, m, d]
x0 = [0; 0; 0; 0];
x_guess = [0; 0; 0; 0; 2];
nx = length(x_guess); % 4 states, 3 paramters
ny = length(g(x_guess)); % x and theta
P0 = diag([0; 0.00001; 0; 0.00001; 0.1]);
u0 = 0;
nu = length(u0);

Q = diag([0; 0.00001; 0; 0.00001; 0.00001]); % Model uncertainty
R = 0.001*eye(ny); % Measurement uncertainty


function dx = cartpend(x,u)
% Adapted from code by Steve Brunton
% x contains state, input and parameters
% x = [x;
%     x_dot;
%     theta;
%     theta_dot;
%     L;
%     m;
%     d;
%     M;]

% Parameters
m = 2;% x(6); % 2 actual value
M = 4; % 4
L = x(5); % 1
g = -9.81;
d = 5; % x(7); % 5

Sx = sin(x(3));
Cx = cos(x(3));
D = m*L*L*(M+m*(1-Cx^2));

D_min = 1e-5;
if D < D_min
    D = D_min;
end

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

function J = jaccsd(f,x,u) % ??? Maybe should use simbolic diff for more exact
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



