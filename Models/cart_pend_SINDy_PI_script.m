%% Performs SINDy-PI on pendulum on a cart system.
% Based on code by Eureka Kaiser
% github.com/dynamicslab/SINDy-PI

%% First try with simple system called "system_ODE"
close all;
% clear all;

tic;

%% Read Data

% Load: x0 and out from simulation
load('rational_toy_with_input_data_1.mat')

% x0 = [0.1; -0.2; 0.3];  % Initial condition
% n = length(x0);  % Number of states
% tspan = [0.01:0.01:10];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
% [t,x] = ode45(@(t,x) rational_toy_ODE(t,x), tspan, x0, options);

% Calculate derivatives
% ??? Change to calculate dx with total variation derivative
% x_dot = rational_toy_ODE(t,x');

t = out.tout;
X = out.x.Data;
X_dot = out.x_dot.Data;
U = out.u.Data;

n = size(X,2); % Number of states

% Add noise to measurements
sigma = 0.001; % Magnitude of noise (max 0.01 so far)
sigma = 0;
X       = X + sigma*randn(size(X));
X_dot   = X_dot + sigma*randn(size(X_dot)); 

% Plot measurements
plot(t,X);
figure

% Theta to compute function for x_dot
Theta_X = Theta(X,U);
num_functions = size(Theta_X,2)*2; % Number of funxtions in Theta(x,xdot) library

%% Find best model for each state
Xi = []; % Stores final Xi, which each column a xi per state
model_errors = []; % error of model for each state
model_lambdas = []; % lambda used for model of each 
model_column_guess = []; % Column used as guess for model of each state

warning('off','MATLAB:rankDeficientMatrix'); % Do not warn about rank deficiency

for i = 1:n % Loop through all states, i
    Theta_i = [Theta_X, diag(X_dot(:,i))*Theta_X]; % Theta used for x1 only
    
    min_metric = Inf; % Store minimum model 2-norm error
    best_lambda = Inf; % Store best lambda
    best_column = Inf; % Store best guess of column
    best_xi = []; % Store best xi for model
   
    lambda_list = logspace(-3,-1,10); % List of lambdas to try
   
    for lambda = lambda_list % Test each lambda in list
        tic_lambda=tic();
        for j = 1:1:size(Theta_i,2)/2+1 % Guess a column in Theta (only for numerator)

            Theta_rm = remove_column(Theta_i,j); % Remove column being guessed
             
            % sparsifyDynamics hogs the time because of backslash/regression:
            xi = sparsifyDynamics(Theta_rm,Theta_i(:,j),lambda); % Sequential LS
            
            % Calculate 2-norm error without test data, penalised for number of parameters
            num_terms = nnz(xi)+1; % Number of non-zero terms in model
            error1 = norm(Theta_i(:,j) - Theta_rm*xi, 1)/norm(Theta_i(:,j),1);
            metric = error1*num_terms^2; % Metric used to compare candidate models
            % Mertric promotes sparsity
            
            % Plot error vs #terms of all candidate models
            subplot(2,n,i), semilogy(num_terms,error1, 'x'), hold on;
            subplot(2,n,i+n), semilogy(num_terms,metric, 'x'), hold on;          

            % ??? Maybe try change to AIC
            
            % Update best_xi if metric is smaller than smallest metric yet
            if metric < min_metric
                best_xi = [xi(1:j-1); -1; xi(j:end)];  % Insert -1 into position j for removed column of Theta
                min_metric = metric; % Update min_error value
                best_error1 = error1; % Save 2-norm error of original data
                best_lambda = lambda; % Store lambda value used for this column
                best_column = j; % Store column used for this guess
            end
            
        end % End: for each column, j, in Theta
        
    end % End: for each lambda in lambda_list
    
    subplot(2,n,i), plot(nnz(best_xi),best_error1, 'o')
    ylabel('data error');
    grid on;
    hold off;
    
    subplot(2,n,i+n), plot(nnz(best_xi),min_metric, 'o')
    ylabel('metric');
    grid on;
    hold off;
    
    % Append model to Xi for this state
    Xi(:,i) = best_xi;
    model_errors(:,i) = min_metric;
    model_lambdas(:,i) = best_lambda;
    model_column_guess(:,i) = best_column;
    
end % End: for each state, i

% Find max/min y scales for scatter plots
subplot(2,n,1), y_scale = ylim;
y_max = max(y_scale);
y_min = min(y_scale);
for i = 2:n*2
    subplot(2,n,i), y_scale = ylim;
    if max(y_scale) > y_max
        y_max = max(y_scale);
    end
    if min(y_scale) < y_min
        y_min = min(y_scale);
    end
end

% Set equal y_scales for scatter plots
y_min = y_min*1e-1; % Give viewing space before limit
y_max = y_max*1e1;
for i = 1:n*2
    subplot(2,n,i), ylim([y_min y_max]);
end

for i = 1:n*2
    subplot(2,n,i), ylim
end

toc

%% Compare real Xi to model Xi with bar plots
figure;
for i = 1:n
    subplot(1,n,i), bar3([Xi(:,i), real_Xi(:,i)]);
end

%% Visualise Xi
x_names = {'x1', 'x2', 'x3', 'u', 'sin(x1)'};
% x_names = [x_names, {'1'}, x_names];
vis_Xi = visualize_Xi(x_names, Xi, 2);

model_errors
model_lambdas
model_column_guess % Guessed column of Theta that worked best


%% Validatation
% Run model on new data and compare to actual measurements

% Load xo, out from a simulation
load('rational_toy_with_input_data_2.mat')
tspan = out.tout;
X = out.x.Data;
U = out.u.Data;

% Generate data with SINDY-PI model
% [t_hat,x_hat] = ode45(@(t,x) SINDy_PI_ODE(t, x, U(), Xi), t_data, x0, options);
x_hat(1,:) = x0;
t_hat(1,:) = 0;

% Solve for small intervals with constant u
for i=2:size(U,1)
    [t_1,x_1] = ode45(@(t_1,y_1) SINDy_PI_ODE(t_1,y_1,U(i-1,:),Xi), tspan(i-1:i,1), x0);
    x_hat(i,:) = x_1(end,:);
    t_hat(i,:) = t_1(end,:);
    x0 = x_hat(i,:)';
end

% PLot simulation data vs model data
figure
plot(tspan,X); hold on;
plot(t_hat,x_hat,'--', 'LineWidth', 1); hold off;

toc; % Display execution time

warning('on','MATLAB:rankDeficientMatrix'); % Switch on warning for other scripts

function dx = rational_toy_ODE(t,x,u)
    % 2 states
    x1 = x(1,:);
    x2 = x(2,:);
    x3 = x(3,:);
    
    dx = [
            - x1 + 3*x3 - x2.^2 + u;
            (10*x1  - x2 - 2*x1.*x3)./(1.1 + sin(x1) + x1.^2);
            (x1.*x2 - 3*x3)
    ];

end

function Theta = Theta(X, U)
%     x1 = X(:,1);
%     x2 = X(:,2);
%     x3 = X(:,3);
%     
%     Theta = [ones(samples,1), x1, x2, x3, x1.*x2, x1.*x3, x1.^2, x2.^2, x3.^2];
    X = [X, U, sin(X(:,1))];
    Theta = [ones(size(X,1),1), X];
    for i = 1:size(X,2)
        for j = i:size(X,2)
            Theta = [Theta, X(:,i).*X(:,j)];
        end
    end
    
end

function x_dot = SINDy_PI_ODE(t, x, u, Xi)
    % ODE created by SINDy_PI model defined by Theta and Xi
    % x is column state vector
    % x_dot is column derivative vector
    n = length(x); % Number of states
    
    X = x'; % Row vector fits SINDY format  
    Theta_X = Theta(X,u); % Calculate library of functions once
    X_dot = zeros(size(X)); % Allocate space once
    for i = 1:n % For each state
        xi = Xi(:,i); 
        xi_num = xi(1:end/2); % Top half of xi for model numerator
        xi_den = xi(end/2+1:end); % Bottom half of xi for model denomenator
        X_dot(:,i) = -(Theta_X*xi_num)./(Theta_X*xi_den);
    end
    
    x_dot = X_dot'; % Back to column vector
end

function Theta_rm = remove_column(Theta,column)
    % Remove the column of Theta given by index: column
    num_columns = size(Theta,2);
    Theta_rm = Theta(:,(1:num_columns)~=column);
end

function Xi = sparsifyDynamics(Theta,dXdt,lambda)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz

    % compute Sparse regression: sequential least squares
    Xi = Theta\dXdt;  % initial guess: Least-squares
    Xi_prev = Xi;
    
    % lambda is our sparsification knob.
    for k=1:10
        small_indexes = (abs(Xi)<lambda);   % find small coefficients
        Xi(small_indexes)=0;                % set small coeffs to 0 (threshold)

        big_indexes = ~small_indexes;
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(big_indexes,:) = Theta(:,big_indexes)\dXdt;
        
        if(Xi == Xi_prev)
            % To save computations, almost by half:
            % If Xi already converged, then exit loop
            break;
        end
        Xi_prev = Xi; % Save previous Xi
    end
end