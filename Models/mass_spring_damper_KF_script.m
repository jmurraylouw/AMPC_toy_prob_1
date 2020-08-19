% Read simulation data
u_data = out.u.Data';
x_data = out.x.Data';
y_data = out.y.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(y_data);

% System definition
m = 1;
b = 0.01;
k = 5;

A= [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = t(2)-t(1);
sys_c = ss(A,B,C,D);
sys_d = c2d(sys_c, Ts);
[F,G,H,D] = ssdata(sys_d);

x_hat_data = zeros(nx, n_time);

sigma_a = 0.1; % Std dev of acceleration/force applied to model
Q = 0.000001*eye(nx); % Model uncertainty
R = 0.01*eye(ny); % Measurement uncertainty

% Initialise
x0 = [0; 0];
P0 = [0, 0; 0, 0];
x_hat_dwork = x0;
P_dwork = P0;

for n = 1:1:n_time-1
    % Measurement
    y = y_data(n);
    u = u_data(:,n); %% ?? Shouldn't u be in y, since it is measured?
    
    % Get saved data
    x_hat = x_hat_dwork;
    P = P_dwork;
    
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Update estimate with measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,n) = x_hat;
    
    x_hat = F*x_hat + G*u; % Extrapolate state
    P = F*P*F' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end
figure
plot(t, x_data(1,:)); hold on;

plot(t, x_hat_data(1,:));
%plot(t, y_data);
hold off;
legend('Actual', 'Estimate', 'Measured')






