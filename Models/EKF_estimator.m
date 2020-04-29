% Read simulation data
u_data = out.f.Data';
x_data = out.x.Data';
y_data = out.y.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(y_data);
nx
nu
ny
n_time

f = @cartpend; % Function handle
x = [0; 0; 0; 0];

% System definition
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



