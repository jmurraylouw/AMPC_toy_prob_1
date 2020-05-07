%  Single mass spring damper
% Try use Last squares to find parameters of mass spring damper
% Based on: Pietersen, System Identification for Fault Tolerant Control of Unmanned Aerial Vehicles

% ??? Why does anylital solution say x(n) is also dependant on u(n)

disp("Least Squares Parameter Estimation")
disp("----------------------------------")

% Original paramters
m = 1;
b = 0.5;
k = 5;

% Read simulation data
try
    u_data = out.u.Data';
    x_data = out.x.Data';
    y_data = out.y.Data';
    t = out.x.Time';
catch
    disp("No data to MSD read. Using saved data.")
    load('mass_spring_damper_Square_wave_input_Measure_position.mat')
end

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(y_data);

% State space system
A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
T = t(2)-t(1);
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c, T);
[F,G,H,D] = ssdata(sys_d);

% Assign memory
N = n_time-2; % Number of measurements
n_p = 4; % Number of Regressors (columns)
z = zeros(N,1); % Measurement vector (z = y + v)
X = zeros(N,n_p); % Regressors
theta = zeros(n_p,1); % Parameter vector

% Populate matrixes
z = y_data(3:end)';
X = [y_data(2:end-1)', y_data(1:end-2)', u_data(2:end-1)', u_data(1:end-2)'];

% Least squares
theta_hat = inv(X'*X)*X'*z

x_hat_data = zeros(n_time,2);
% Initial condition
x_hat_data(1,1) = 0; 
x_hat_data(1,2) = 0;

for n = 3:1:n_time
    x = [x_hat_data(1,n-1), x_hat_data(1,n-2), u_data(n-1), u_data(n-2)];
    x_hat_data(1,n) = x*theta_hat;
end

figure
plot(t,x_data(1,:))
hold on;
plot(t,x_hat_data(1,:), '--')
plot(t,x_data(1,:)-x_hat_data(1,:))
hold off
legend("Actual x", "Data Model x", "Error")
title("Actual data vs data-driven model")

% Mean squared Model error
MSE_data = mean((x_data(1,:)-x_hat_data(1,:)).^2)

% Compare to analtyical discrete model
% Tuskin model     
theta = [...
     ((8*m) - (2*T^2*k))/(k*T^2 + 2*b*T + 4*m);
     ((2*T*b)/(k*T^2 + 2*b*T + 4*m) - (4*m)/(k*T^2 + 2*b*T + 4*m) - (T^2*k)/(k*T^2 + 2*b*T + 4*m));
     
     ((2*T^2)/(k*T^2 + 2*b*T + 4*m));
     (2*T^2/(k*T^2 + 2*b*T + 4*m))]
     
x_hat_data = zeros(n_time,2);
% Initial condition
x_hat_data(1,1) = 0; 
x_hat_data(1,2) = 0;

for n = 3:1:n_time
    x = [x_hat_data(1,n-1), x_hat_data(1,n-2), u_data(n-1), u_data(n-2)];
    x_hat_data(1,n) = x*theta;
end

figure
plot(t,x_data(1,:))
hold on;
plot(t,x_hat_data(1,:), '--')
plot(t,x_data(1,:)-x_hat_data(1,:))
hold off
legend("Actual x", "Analytic Model x", "Error")
title("Actual data vs analytical model")

% Mean squared Model error
MSE_analytic = mean((x_data(1,:)-x_hat_data(1,:)).^2)

% Solve for parameters
syms m k b
eqn_theta = [...
     ((8*m) - (2*T^2*k))/(k*T^2 + 2*b*T + 4*m);
     ((2*T*b - (4*m) - (T^2*k))/(k*T^2 + 2*b*T + 4*m));
     
     ((2*T^2)/(k*T^2 + 2*b*T + 4*m));
     (2*T^2/(k*T^2 + 2*b*T + 4*m))] ...
     == theta_hat;

[num_L, denom_L] = numden(lhs(eqn_theta)); % extract left numerator and denomenator
[num_R, denom_R] = numden(rhs(eqn_theta));
% Flatten all fractions
eqn_theta = 0 == num_L.*denom_R - num_R.*denom_L;

Coef = zeros(n_p, 4);

old = [1 m b k];
new = [1 2 3 4]; % Index where each symbol should be
term_indexes = subs(term, old, new);
for row = 1:1:n_p
    [coef,term] = coeffs(rhs(eqn_theta(row)));
    term_indexes = subs(term, old, new);
    for index = 1:1:length(term_indexes)
        col = term_indexes(index);
        Coef(row,col) =  coef(index);
    end
end

% Solve overdetermined equations to get parameters
% A*parameters = B
A = Coef(:,2:end);
B = -Coef(:,1);
parameters = A\B; % Least squares solution to overdetermined system
m_hat = parameters(1)
b_hat = parameters(2)
k_hat = parameters(3)

m = 1.1;
b = 0.2;
k = 3;


m = 1;
b = 0.5;
k = 5;






