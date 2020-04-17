% Basesed on algorithm explained by cs.tut.fi/~tabus/course/ASP/LectureNew10.pdf
% Code written by: JM Louw
disp("START")
disp("-----")
f       = out.f.data; % input data/force values
x       = out.x.data; % output data/distance values
N       = max(size(x)); % Number of time steps
eps     = zeros(1,N); % Error between predicted and actual output
y       = zeros(1,N); % Predicted output

M = 6; % Number of parameters in Theta

% Initialise
n0      = 3; % Time step to start algorithm at, because of (k-2)
Theta   = zeros(M,1);     % Inititialise parameter vectors
P       = 100*eye(M);   % Initialise large P
lambda  = 0.99;         % forgetting factor

% No need to remeber previous k or P,
% only update k and P for next use.

% RLS
for n = n0:1:N
    phi     = [f(n); f(n-1); x(n-1); x(n-2); eps(n-1); eps(n-2)]; % 2nd order model, therefore 2 noise parameters 
    pi      = phi'*P;
    gamma   = lambda + pi*phi;
    k       = pi'/gamma;
    y(n)    = Theta'*phi;
    eps(n)  = x(n) - y(n);
    Theta   = Theta + k*eps(n); % Extimate next paramter vector

    P_prime = k*pi;
    P       = 1/lambda*(P - P_prime); % Update P for next use
end

error_ave = sum(abs(eps))/2;
figure;
plot(y); hold on
plot(x)
plot(eps, 'k'); hold off

legend('x predicted (y)', 'x actual', 'error (eps)')
title(strcat('l = ', num2str(lambda), '; noise ps = ', int2str(M-4), '; err ave = ', num2str(error_ave)))
disp("END")
disp("---")