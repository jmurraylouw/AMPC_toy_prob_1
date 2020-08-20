% Setup workspace for EKF of mass spring damper

% Dimensions
nx = 2;
nu = 1;
ny = 1;

% Linear msd system definition
m = 1;
b = 0.01;
k = 5;

A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;

% Sample time of estimator
Ts = 0.01;

% Initialise
x0 = [0; 0];
P0 = [0, 0; 0, 0];
u0 = 0;

Q = 0.00001*eye(nx); % Model uncertainty
R = 0.001*eye(ny); % Measurement uncertainty

% Function handles
f = @nl_msd; % System function handle
g = @measure_msd; % Measurement function handle

function dx = msd(x,u)
    % LINEAR MASS SPRING DAMPER
    % INPUT:
    % x = state vector [x, x_dot]
    % u = input vector [f]
    % 
    % OUTPUT:
    % dx = derivative of state vector [x_dot, x_dotdot]

    m = 1; % Mass
    b = 0.01; % Damper coef
    k = 5; % Coef of spring

    dx = zeros(2,1); % Assign memory
    dx(1,1) = x(2);
    dx(2,1) = 1/m*(-k*x(1) - b*x(2) + u);
end

function dx = nl_msd(x,u)
    % NON-LINEAR MASS SPRING DAMPER
    % INPUT:
    % x = state vector [x, x_dot]
    % u = input vector [f]
    % 
    % OUTPUT:
    % dx = derivative of state vector [x_dot, x_dotdot]

    m = 1; % Mass
    b = 0.01; % Damper coef
    k = 5; % Coef of non-linear spring

    dx = zeros(2,1); % Assign memory
    dx(1,1) = x(2);
    dx(2,1) = 1/m*(-k*x(1)^3 - b*x(2) + u);
end

function y = measure_msd(x,u)
    % Measurement function
    % g = @measure_msd
    % C = Jacobian of g
    y = x(1);
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
    h = n*eps;
    for k=1:n
        x1 = x;
        x1(k) = x1(k)+ h*1i;
        J(:,k) = imag(f(x1,u))/h;
    end
end



