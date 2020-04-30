% Read simulation data
x_data = out.x.Data';
y_data = out.y.Data';
u_data = out.u.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[ny, n_time] = size(y_data);
[nu, n_time] = size(u_data);

% System definition
m = 1;
b = 0.01;
k = 5;

A= [0, 1; -k/m, -b/m]
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = t(2)-t(1);

% Initialise
x0 = [0; 0];
P0 = [0, 0; 0, 0];
x_hat = x0;
P = P0;

Q = 0.000001*eye(nx); % Model uncertainty
R = 0.01*eye(ny); % Measurement uncertainty

% Extrapolate
f = @msd; % System function handle
g = @measure_msd; % Measurement function handle
H = jaccsd(g,x_hat,0); %H = C if measurement system is linear
% ??? Change this if f is dependent on time.
x_hat_dwork = x_hat + f(x_hat,0)*Ts; % Numeric integration to extrapolate state

F = jaccsd(f,x_hat,0); % Calculate Jacobian of continuous system
Phi = eye(nx) + F*Ts + 1/2*(F*Ts)^2; % ??? where is this from? 2nd order Taylor expansion? (continuous to discrete)
P_dwork = Phi*P*Phi' + Q; % Extrapolate uncertainty

x_hat_data = zeros(nx, n_time); % Assign memory beforehand

% Apply EKF at every timestep
for n = 1:1:n_time-1
    % Measurement
    y = y_data(n);
    u = u_data(n);
    
    % Get saved data
    x_hat = x_hat_dwork;
    x_hat_dwork;
    P = P_dwork;
    
    % Update
    H = jaccsd(g,x_hat,0); % Linearise measurement function
    K = (P*H')/(H*P*H' + R); % Compute Kalman gain (b*inv(A) -> b/A)
    x_hat = x_hat + K*(y - H*x_hat); % Update estimate with measurement
    KH_term = (eye(nx) - K*H);
    P = KH_term*P*KH_term' + K*R*K'; % Update estimate uncertainty
    
    % Output
    x_hat_data(:,n) = x_hat;
    
    % Extrapolate for next time step
    x_hat = x_hat + f(x_hat,u)*Ts; % Numeric integration (extrapolate state)

    F = jaccsd(f,x_hat,0); % Calculate Jacobian of continuous system
    Phi = eye(nx) + F*Ts + 0.5*(F*Ts)^2; % 2nd order Taylor expansion (continuous to discrete)
    P = Phi*P*Phi' + Q; % Extrapolate uncertainty
    
    % Save to Dwork
    x_hat_dwork = x_hat;
    P_dwork = P;
    
end

figure
plot(t,x_data(1,:))
hold on
plot(t,x_hat_data(1,:))
hold off

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

function dx = msd(x,u)
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
    dx(2,1) = 1/m*(-k*x(1) - b*x(2) + u);
end

function y = measure_msd(x,u)
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



