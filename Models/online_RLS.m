% Basesed on algorithm explained by cs.tut.fi/~tabus/course/ASP/LectureNew10.pdf
% Code written by: JM Louw
disp("START")
disp("-----")
f_data       = out.f.data; % input data/force values
x_data       = out.x.data; % output data/distance values
N            = max(size(x_data)); % Number of time steps
eps_data     = zeros(1,N); % Error between predicted and actual output
y_data       = zeros(1,N); % Predicted output

M = 6; % Number of parameters in Theta

% Initialise
lambda  = 0.99;         % forgetting factor

% No need to remeber previous k or P,
% only update k and P for next use.

% RLS
prev_data = zeros(M-1,1);
Theta   = zeros(M,1); % Inititialise parameter vectors, Time steps downwards
P       = 100*eye(M); % Initialise large P
P2      = reshape(P, M*M, 1); % Form a vector to store in Dwork

% Assign to memory/Dwork
prev_data_dwork = prev_data;
Theta_dwork     = Theta;
P2_dwork        = P2;
    
for n = 1:1:N
    disp(n)
    disp("-----")
    f       = f_data(n);
    x       = x_data(n);

    prev_data = prev_data_dwork;
    Theta     = Theta_dwork;
    P2        = P2_dwork;

    P       = reshape(P2, M, M); % Reshape vector P2 in matrix P
    % prev_data = [f(n-1); x(n-1); x(n-2); eps(n-1); eps(n-2)]
    phi     = [f; prev_data]; % 2nd order model, therefore 2 noise parameters 
    pi      = phi'*P; % Use M previous values of U as well
    gamma   = lambda + pi*phi;
    k       = pi'/gamma;
    y       = Theta'*phi;
    eps     = x - y;
    Theta   = Theta + k*eps; % Extimate next w vector

    P_prime = k*pi;
    P       = 1/lambda*(P - P_prime); % Update P for next use
    P2          = reshape(P, M*M, 1); % Form a vector to store in Dwork   
    prev_data   = [f; x; prev_data(2); eps; prev_data(4)];
    
    y_data(n) = y;
    eps_data(n) = eps;
    
    prev_data_dwork = prev_data;
    Theta_dwork = Theta;
    P2_dwork = P2;
    P       = 1/lambda*(P - P_prime); % Update P for next use
end

error_ave = sum(abs(eps_data))/2;
figure;
plot(y_data); hold on
plot(x_data)
plot(eps_data, 'k'); hold off

legend('x predicted (y)', 'x actual', 'error (eps)')
title(strcat('l = ', num2str(lambda), '; noise ps = ', int2str(M-4), '; err ave = ', num2str(error_ave)))
disp("END")
disp("---")