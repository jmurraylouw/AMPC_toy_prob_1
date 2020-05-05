%  Single mass spring damper
% Try use Last squares to find parameters of mass spring damper
% Based on: Pietersen, System Identification for Fault Tolerant Control of Unmanned Aerial Vehicles

m = 1;
b = 0.5;
k = 5;

%%
% Read simulation data
u_data = out.u.Data';
x_data = out.x.Data';
y_data = out.y.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(y_data);

% State space system
A= [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = t(2)-t(1);
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c, Ts);
[F,G,H,D] = ssdata(sys_d);

% Assign memory
N = n_time-2; % Number of measurements
n_p = 4; % Number of parameters
z = zeros(N,1); % Measurement vector (z = y + v)
X = zeros(N,n_p); % Regressors
theta = zeros(n_p,1); % Parameter vector

% Populate matrixes
z = y_data(3:end)';
X = [ones(N,1), y_data(2:end-1)', y_data(1:end-2)', u_data(2:end-1)'];


% Least squares
theta_hat = inv(X'*X)*X'*z

x_hat_data = zeros(n_time,2);
% Initial condition
x_hat_data(1,1) = 0; 
x_hat_data(1,2) = 0;

for n = 3:1:n_time
    x = [1, x_hat_data(1,n-1), x_hat_data(1,n-2), u_data(n-1)];
    x_hat_data(1,n) = x*theta_hat;
end

plot(t,y_data)
hold on;
plot(t,x_hat_data(1,:), '--')
hold off








