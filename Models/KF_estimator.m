u_data = out.f.Data';
x_data = out.x.Data';
z_data = out.z.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(z_data);

x_hat_data = zeros(nx, n_time); % a priori prediction

Q = zeros(nx, nx);
R = 0.01;

% System definition
m = 1;
b = 0.5;
k = 5;

A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = t(2)-t(1);
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c, Ts)
[A,B,C,D] = ssdata(sys_d)

% Initialise
x_hat = [0; 0];
P = [0, 0; 0, 0];

x_hat_dwork = A*x_hat + B*u % Extrapolate state
P_dwork = A*P*A' + Q; % Extrapolate uncertainty

 
for k = 1:1:n_time-1
    % Measurement
    z = z_data(k);
    u = u_data(:,k);
    
    % Get saved data
    x_hat = x_hat_dwork;
    P = P_dwork;
    
    K = P*C'*inv(C*P*C' + R); % Compute Kalman gain
    x_hat = x_hat + K*(z - C*x_hat); % Update estimate with measurement
    P = (1 - K*C)*P*(1 - K*C)' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,k) = x_hat;
    
    x_hat = A*x_hat + B*u; % Extrapolate state
    P = A*P*A' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end

plot(t, x_data(1,:)); hold on;
plot(t, z_data);
plot(t, x_hat_data(1,:)); hold off;
legend('Actual', 'Measured', 'Estimate')






