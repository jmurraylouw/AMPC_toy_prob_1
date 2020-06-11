% mass_spring_damper_partial_D_SINDy_script

%% Gather data
%% Generate Data
x0 = [-1; 0.1];  % Initial condition
x0_valid = [2; 2];

nx = length(x0); % Number of states

tspan = 0:0.01:20;
tspan_valid = [0:0.01:40];

% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
ode = @toy_ODE; % ODE system to use
[t,x] = ode45(@(t,x) ode(t,x), tspan, x0); % Trainging data
[t_valid,x_valid] = ode45(@(t,x) ode(t,x), tspan_valid, x0_valid); % Validation data

sigma  = 0.1; % Standard deviation of noise 
x = x + sigma*rand(size(x)); % Add noise
Y = x(:,1)'; % only measure position

% figure(1), plot(t,x);
% pause

%% Step 1: Construct snapshot matrix
num_delays = 50;
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
plot(V_tilde);
pause

%% Step 3: SVD and truncate output matrix, X2
[U,S,V] = svd(X2,'econ');

% semilogy(diag(S), 'o')
r = 2; % truncate order
% Reduced rank of Omega svd
U_hat = U(:, 1:r); % Truncate SVD matrixes of Omega
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);
% plot(V_hat)

%% Step 4: SINDy regression
lambda = 0.001; % lambda is our sparsification knob.
polyorder = 1;
Theta_V_tilde = Theta(V_tilde,polyorder);
Xi = [];
for i = 1:r
    Xi(:,i) = sparsify_dynamics(Theta_V_tilde,V_hat(:,i),lambda)
end
V_est = Theta_V_tilde*Xi;
plot(V_hat); hold on
plot(V_est, '--', 'LineWidth', 2)

%% Step 5: Run model
% 5.0: Setup Hankel with initial condition data
x_hat = zeros(r,length(t));
H = [];
for i = 1:num_delays
    H = [Y(:,i:window-1+i); H];
end

% 5.0: Setup Hankel with initial condition data

for k = num_delays:length(t)
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








