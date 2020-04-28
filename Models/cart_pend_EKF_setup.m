% sfunc_EKF(A,B,C,D,Ts,Q,R,x0,P0)

Ts = 0.01;

nx = 4; % [x; x_dot; theta; theta_dot]
nu = 1; % f (horizontal force on cart)
ny = 4;

m = 1;
M = 5;
L = 2;
g = -10;
d = 1;

s = -1; % pendulum up (s = 1), pend down (s = -1)

% System matrixes
A = [0 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0];

B = [0; 1/M; 0; s*1/(M*L)];

C = eye(4);
%C = [1 0 0 0;
%     0 0 1 0];
D = 0;

% Uncertainties
Q = zeros(nx, nx);
R = 0.01*eye(ny);

% Initial estimates
x0 = [0; 0; 0; 0];
P0 = zeros(nx, nx);

disp('cart_pend KF data loaded')