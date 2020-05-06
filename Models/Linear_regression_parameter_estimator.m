%  Single mass spring damper
% Try use Last squares to find parameters of mass spring damper
% Based on: Pietersen, System Identification for Fault Tolerant Control of Unmanned Aerial Vehicles

m = 1;
b = 0.5;
k = 5;


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
T = t(2)-t(1);
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

figure
plot(t,y_data)
hold on;
plot(t,x_hat_data(1,:), '--')
plot(t,y_data-x_hat_data(1,:))
hold off
legend("Actual x", "Data Model x", "Error")
title("Actual data vs data-driven model")

% Mean squared Model error
MSE_data = mean((y_data-x_hat_data(1,:)).^2)

% Compare to analtyical discrete model
% 1/s^2 = Tz/(z-1)^2
theta = [0;
         (2*m+b-k*Ts)/(m+b);
         -m/(m+b);
         Ts/(m+b)]; % From analytical derivation. Differential eq -> Laplace -> Z-transform -> Difference eq

% Tuskin model     
theta = [0;
     ((8*m) - (2*T^2*k))/(k*T^2 + 2*b*T + 4*m);
     ((2*T*b)/(k*T^2 + 2*b*T + 4*m) - (4*m)/(k*T^2 + 2*b*T + 4*m) - (T^2*k)/(k*T^2 + 2*b*T + 4*m));
     (T^2/(k*T^2 + 2*b*T + 4*m));
     ((2*T^2)/(k*T^2 + 2*b*T + 4*m));
     (T^2/(k*T^2 + 2*b*T + 4*m))]

     
x_hat_data = zeros(n_time,2);
% Initial condition
x_hat_data(1,1) = 0; 
x_hat_data(1,2) = 0;

for n = 3:1:n_time
    x = [1, x_hat_data(1,n-1), x_hat_data(1,n-2), u_data(n), u_data(n-1), u_data(n-2)];
    x_hat_data(1,n) = x*theta;
end

figure
plot(t,y_data)
hold on;
plot(t,x_hat_data(1,:), '--')
plot(t,y_data-x_hat_data(1,:))
hold off
legend("Actual x", "Analytic Model x", "Error")
title("Actual data vs analytical model")

% Mean squared Model error
MSE_analytic = mean((y_data-x_hat_data(1,:)).^2)





