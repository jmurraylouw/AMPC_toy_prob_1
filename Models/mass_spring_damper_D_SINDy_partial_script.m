% mass_spring_damper_partial_D_SINDy_script
clear all;

%% Gather data
%% Generate Data
x0 = [-1; 0.1];  % Initial condition
x0_valid = [2; 2];

nx = length(x0); % Number of states

tspan = 0:0.01:20;
tspan_valid = 0:0.01:10;

% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
ode = @toy_ODE; % ODE system to use
[t,x] = ode45(@(t,x) ode(t,x), tspan, x0); % Trainging data
[t_valid,x_valid] = ode45(@(t,x) ode(t,x), tspan_valid, x0_valid); % Validation data

sigma  = 0; % Standard deviation of noise 
x = x + sigma*rand(size(x)); % Add noise

ny = 1; % Number of measurements
Y = x(:,1)'; % only measure position
Y_valid = x_valid(:,1)'; % only measure position, for validation data

% figure(1), plot(t,x);
% pause

%% Step 1: Construct snapshot matrix
num_delays = 4;
window = size(x,1)-num_delays-1;
H = [];
for i = 1:num_delays
    H = [Y(:,i:window-1+i); H];
end

na = size(H,1); % Length of augmented state vector

X2 = H(:, 2:end); % Y advanced 1 step into future
X = H(:, 1:end-1); % Cut off last sample

%% Step 2: SVD and truncate input matrix, X
[U,S,V] = svd(X,'econ');
% semilogy(diag(S), 'o')
p = 2; % truncate order
% Reduced rank of Omega svd
U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
S_tilde = S(1:p, 1:p);
V_tilde = V(:, 1:p);
% plot(V_tilde); hold on
% pause

%% Step 3: SVD and truncate output matrix, X2
[U,S,V] = svd(X2,'econ');

% semilogy(diag(S), 'o')
r = 2; % truncate order
% Reduced rank of Omega svd
U_hat = U(:, 1:r); % Truncate SVD matrixes of Omega
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);
% plot(V_hat, 'k--'); hold off
% pause

%% Step 4: SINDy regression
lambda = 0.001; % lambda is our sparsification knob.
polyorder = 1;
Theta_V_tilde = Theta(V_tilde,polyorder);
Xi = [];
for i = 1:r
    Xi(:,i) = sparsify_dynamics(Theta_V_tilde,V_hat(:,i),lambda);
end
V_est = Theta_V_tilde*Xi;
% plot(V_hat); hold on
% plot(V_est, '--', 'LineWidth', 2)

%% Step 5: Run model
tic;
% 5.0: Setup initial condition data
start_k = 200; % Minimum k in x(k) from which x(k+1) is predicted
assert(num_delays*2 <= start_k) % Must be at least a square matrix
width_H0 = start_k - num_delays + 1; % Number of columns in initial Hankel

x_hat = zeros(ny,length(t));
x_hat(:,1:start_k) = Y_valid(1,1:start_k); % Initial known data

H = [];
for i = 1:num_delays
    H = [x_hat(:,i:width_H0-2+i); H];
end

for k = start_k:length(t_valid)-1

    % 5.1: Setup Hankel with previous data
    new_column = [x_hat(:,k); H(1:end-1,end)]; % Add new data entry to Hankel column
    H = [H, new_column]; % Add column to Hankel matrix

    % 5.2: SVD of Hankel matrix
    [U,S,V] = svd(H);
    U_tilde = U(:, 1:p); % Truncate SVD matrixes 
    S_tilde = S(1:p, 1:p);
    V_tilde = V(:, 1:p);

    % 5.3: Predict next x
    Theta_V_tilde = Theta(V_tilde,polyorder);
    V_hat = Theta_V_tilde*Xi;
%     plot(V_tilde); hold on
%     pause
    
%     U_hat = U_tilde; % Approximate
%     S_hat = S_tilde; % Approximate
    X2_hat = U_hat*S_hat*V_hat';
    x_hat(k+1) = X2_hat(1,end); % predicted x
    
end

plot(x_valid(:,1)); hold on
plot(x_hat, 'k--'); 
plot(start_k, x_hat(start_k), 'ro')
hold off

toc   
%% Local functions

function dx = toy_ODE(t,x)
    m=1;
    b=0.1;
    k=15;
    dx = [
            x(2);
            -b/m*x(2) - k/m*x(1); % Linear spring
    ];
end

function Theta_X = Theta(X, polyorder)
     
%     Poly order =  Highest order of polynomial term in library
    X = [X];
    
    n = size(X,2); % number of pseudo states
    
    % Polynomial order 1:
    Theta_X = [ones(size(X,1),1), X];
    
    % Polynomial order 2:
    if polyorder >= 2
        for i = 1:n
            for j = i:n
                Theta_X = [Theta_X, X(:,i).*X(:,j)];
            end
        end
    end
    
    % Polynomial order 3:
    if polyorder >= 3
        for i=1:n
            for j=i:n
                for k=j:n
                    Theta_X = [Theta_X, X(:,i).*X(:,j).*X(:,k)];
                end
            end
        end
    end
    
end

function xi = sparsify_dynamics(Theta_X,x2,lambda)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz

    % compute Sparse regression: sequential least squares
    xi = Theta_X\x2;  % initial guess: Least-squares
    Xi_prev = xi;
    
    % lambda is our sparsification knob.
    k_conv = 20; % k at which Xi converges
    for k=1:k_conv
        small_indexes = (abs(xi)<lambda);   % find small coefficients
        xi(small_indexes)=0;                % set small coeffs to 0 (threshold)

        big_indexes = ~small_indexes;
        % Regress dynamics onto remaining terms to find sparse Xi
        xi(big_indexes,:) = Theta_X(:,big_indexes)\x2;
        
        if(xi == Xi_prev)
            % To save computations, almost by half:
            % If Xi already converged, then exit loop
            k_conv = k;
            break;
        end
        Xi_prev = xi; % Save previous Xi
    end
end








