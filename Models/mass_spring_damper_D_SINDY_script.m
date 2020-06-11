%% Discrete SINDY - Moving-Window of non-linear mass spring damper
% Estimate non-linear system model with moving window 
% of data from non-linear mass spring damper
% Full-state feedback

% Adaptation of code by Steve Brunton

clear all;

%% Generate Data
x0 = [-1; 0.1];  % Initial condition
n = length(x0);  % Number of states
tspan = 0:0.01:20;

% Trainging data
% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
ode = @toy_ODE;
[t,x] = ode45(@(t,x) ode(t,x), tspan, x0);

% Add noise

lambda = 0.001; % lambda is our sparsification knob.
sigma  = 0.01;
x = x + sigma*rand(size(x));
figure(1), plot(t,x);
pause
% Validation data
x0_valid = [2; 2];
tspan_valid = [0:0.01:40];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t_valid,x_valid] = ode45(@(t,x) ode(t,x), tspan_valid, x0_valid);

% Run for different number of delays
num_delays_list = 0:1:5;
MSE = [];
for num_delays = num_delays_list
    start = 20;%num_delays+1; % Index at which data set starts, to give space for delay coordinates
    window = length(t)*0.8; % Size of data sample window used in regression

    % No access to velocity
    X2 = x((start+1:start+window+1),:);  % x(k+1), One timestep ahead

    % Build library and compute sparse regression
    polyorder = 1;
    x_extended = []; % Add other terms to use in Theta

    for i = 0:num_delays
        x_extended = [x_extended, x((start-i:start+window-i),:)];
    end
    
    Xi = [];
    Theta_X = poolData(x_extended,polyorder);  % up to third order polynomials
    for i = 1:n
        Xi(:,i) = sparsifyDynamics(Theta_X,X2(:,i),lambda);
    end
    
    % xi = lasso(Theta_X,X2,'Lambda',lambda);

    x_names = {'k', '(k-1)', '(k-2)', '(k-3)'};
    x_names = {'k'};

    vis_Xi = visualize_Xi(x_names, Xi, polyorder)
    num_delays
    
    % Run model on validation data
    x_hat = zeros(size(x_valid));

    for i = 1:num_delays+1 
        x_hat(i) = x_valid(i); % Set initial conditions
    end

    for k = num_delays+1:length(t_valid)-1
        x_extended = [];
        for i = 0:num_delays
            x_extended = [x_extended, x_hat(k-i,:)];
        end
        Theta_X = poolData(x_extended,polyorder);
        x_hat(k+1,:) = Theta_X*Xi;
    end

    subplot(1,2,1), plot(t_valid,x_valid); hold on
    plot(t_valid,x_hat,'k--','LineWidth', 1); hold off

    % Mean Squared Error on validation data
    subplot(1,2,2), bar3(Xi./sum(abs(Xi)))
    MSE = [MSE, mean(mean(x_valid.^2 - x_hat.^2))];
    MSE(end)
pause

end
figure
semilogy(num_delays_list,MSE)

%%

function dx = toy_ODE(t,x)
    m=1;
    b=0.1;
    k=15;
    dx = [
            x(2);
            -b/m*x(2) - k/m*x(1); % Linear spring
    ];
end

function xi = sparsifyDynamics(Theta_X,x2,lambda)
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


