% Basesed on algorithm explained by cs.tut.fi/~tabus/course/ASP/LectureNew10.pdf
% Code written by: JM Louw
disp("START")
disp("-----")
u_data  = out.u.data'; % input data/force values
x_data	= out.x.data'; % output data/distance values
y_data  = out.y.data'; % output data/distance values
t       = out.tout';    % Output time series

N            = max(size(x_data)); % Number of time steps
eps_data     = zeros(1,N); % Error between predicted and actual output
x_hat_data   = zeros(1,N); % Predicted output

M = 5; % Number of parameters in Theta

% Initialise
lambda  = 0.999;         % forgetting factor

% No need to remeber previous k or P,
% only update k and P for next use.

% RLS
prev_data = zeros(M-1,1);
Theta   = zeros(M,1); % Inititialise parameter vectors
P       = 100*eye(M); % Initialise large P
P2      = reshape(P, M*M, 1); % Form a vector to store in Dwork
Theta_data = zeros(M,N);

% Assign to memory/Dwork
prev_data_dwork = prev_data;
Theta_dwork     = Theta;
P2_dwork        = P2;
    
for n = 1:1:N
    u       = u_data(n); % Input force
    x       = y_data(n); % Measure

    prev_data = prev_data_dwork;
    Theta     = Theta_dwork;
    P2        = P2_dwork;

    P       = reshape(P2, M, M); % Reshape vector P2 in matrix P
    
    % phi = [u(n); u(n-1); u(n-2); x(n-1); x(n-2)]
    % Theta = row vector of parameters
    % x(n) = Theta*phi    
    
    phi     = [u; prev_data]; % 2nd order model, therefore 2 noise parameters 
    pi      = phi'*P; % phi(n)'*P(n-1)
    gamma   = lambda + pi*phi; % lambda + phi(n)'*P(n-1)*phi(n)
    K       = pi'/gamma; % phi(n)'*P(n-1) * inv(lambda + phi(n)'*P(n-1)*phi(n))
    x_hat   = Theta'*phi;
    eps     = x - x_hat;
    Theta   = Theta + K*eps; % Estimate next parameter vector

    P_prime = K*pi;
    P       = 1/lambda*(P - P_prime); % P(n), Update P for next use
    P2      = reshape(P, M*M, 1); % Form a vector to store in Dwork   
    
    % Currently: prev_data = [u(n-1); u(n-2); x(n-1); x(n-2)]
    % Needs to be for next step: prev_data = [u(n); u(n-1); x(n); x(n-1)]
    prev_data   = [u; prev_data(1); x; prev_data(3)]; % Will be used for phi
    
    eps_data(n) = eps;
    Theta_data(:,n) = Theta;
    
    prev_data_dwork = prev_data;
    Theta_dwork     = Theta;
    P2_dwork        = P2;
end


% Plot model
x_hat_data = zeros(1,N)
x_hat_data(1) = 0; % Initial condition

for n = 3:1:N
    phi = [u_data(n); u_data(n-1); u_data(n-2); x_hat_data(n-1); x_hat_data(n-2)];
    x_hat_data(n) = Theta'*phi;
end

eps_data = x_data(1,:) - x_hat_data;

figure;
plot(t,x_hat_data); hold on
plot(t,x_data(1,:), '--')
plot(t,eps_data, 'k'); 
% plot(t,y_data)
hold off

legend('x_hat', 'x actual', 'error (eps)', 'measured')
title(strcat('l = ', num2str(lambda), '; noise ps = ', int2str(M-4), '; err ave = ', num2str(MSE)))

Theta
MSE = mean((eps_data).^2) % Mean Squared Error


disp("END")
disp("---")