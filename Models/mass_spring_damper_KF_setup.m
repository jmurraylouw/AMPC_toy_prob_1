nx = 2;
nu = 1;
ny = 1;

% System definition
m = 1;
b = 0.01;
k = 5;

F= [0, 1; -k/m, -b/m];
G = [0; 1/m];
H = [1, 0];
D = 0;
Ts = 0.01;
sys_c = ss(F,G,H,D);
sys_d = c2d(sys_c, Ts);
[F,G,H,D] = ssdata(sys_d);

x_hat_data = zeros(nx, n_time);

sigma_a = 0.1; % Std dev of acceleration/force applied to model
Q = G*(sigma_a^2)*G'
R = 0.01;

% Initialise
x0 = [0; 0];
P0 = [0, 0; 0, 0];






