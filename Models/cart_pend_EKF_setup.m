% sfunc_EKF(F,G,H,D,Ts,Q,R,x0,P0)

nx = 4; % [x; x_dot; theta; theta_dot]
nu = 1; % f (horizontal force on cart)
ny = 4;

m = 1;
M = 5;
L = 2;
g = -10;
d = 10;

s = -1; % pendulum up (s = 1), pend down (s = -1)

% Continuous System matrixes
A = [0 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0];

B = [0; 1/M; 0; s*1/(M*L)];

C = eye(4);
%C = [1 0 0 0;
%     0 0 1 0];
D = 0;

% Discritize System
Ts = 0.01;
sys_c = ss(A,B,C,D); % Continuous system
sys_d = c2d(sys_c, Ts); % Discrete system
[F,G,H,D] = ssdata(sys_d);

% Uncertainties
sigma_a = 0.1;
Q = G*sigma_a^2*G'; %zeros(nx, nx);
R = 0.00001*eye(ny);

% Initial estimates
x0 = [0; 0; 0.05; 0];
P0 = 0.1*ones(nx, nx);

disp('cart_pend KF data loaded')