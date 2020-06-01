%% Performs SINDy-PI on pendulum on a cart system.
% Based on code by Eureka Kaiser
% github.com/dynamicslab/SINDy-PI

%% First try with simple system called "system_ODE"
tic;
%% Generate Data
clear all;
close all;
x0 = [0.1; -0.2; 0.3];  % Initial condition
x0_test = [-0.5; -0.1; 1]; % Initial condition for test data

n = length(x0);  % Number of states
tspan = [0.01:0.001:10];

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x] = ode45(@(t,x) system_ODE(t,x), tspan, x0, options);
[t_test,x_test] = ode45(@(t,x) system_ODE(t,x), tspan, x0_test, options);

% Calculate derivatives
% ??? Change to calculate dx with total variation derivative
dx = system_ODE(t,x');
dx_test = system_ODE(t_test,x_test');

% Add noise to measurements
sigma = 0.0001; % Magnitude of noise
% sigma =0;
x       = x + sigma*randn(size(x));
x_test  = x_test + sigma*randn(size(x_test));
dx      = dx + sigma*randn(size(dx)); 
dx_test = dx_test + sigma*randn(size(dx_test));

% Plot measurements
subplot(1,2,1), plot(t,x);
subplot(1,2,2), plot(t_test,x_test);
hold off;

% State and derivative matrixes
X = x;
X_dot = dx';

X_test = x_test;
X_dot_test = dx_test';

% Theta to compute function for x_dot
Theta_X = Theta(X);
Theta_X_test = Theta(X_test);
num_functions = size(Theta_X,2)*2; % Number of funxtions in Theta(x,xdot) library

%% Find best model for each state
Xi = []; % Stores final Xi, which each column a xi per state
model_errors = []; % error of model for each state
model_lambdas = []; % lambda used for model of each 
model_column_guess = []; % Column used as guess for model of each state

Xi_r2 = []; % Stores final Xi, which each column a xi per state
model_errors_r2 = []; % error of model for each state
model_lambdas_r2 = []; % lambda used for model of each 
model_column_guess_r2 = []; % Column used as guess for model of each state


warning('off','MATLAB:rankDeficientMatrix'); % Do not warn about rank deficiency

candidate_xi_cell = cell(1,3); % Each element contains matrix with columns as candidate xis for that state
for i = 1:n % Loop through all states, i
    Theta_i = [Theta_X, diag(X_dot(:,i))*Theta_X]; % Theta used for x1 only
    Theta_i_test = [Theta_X_test, diag(X_dot_test(:,i))*Theta_X_test]; % test Theta used for x1 only
    
    min_error = Inf; % Store minimum model 2-norm error
    best_lambda = Inf; % Store best lambda
    best_column = Inf; % Store best guess of column
    best_xi = []; % Store best xi for model

    max_R2adj = 0; % Store minimum model 2-norm error
    best_lambda_r2 = Inf; % Store best lambda
    best_column_r2 = Inf; % Store best guess of column
    best_xi_r2 = []; % Store best xi for model
    
    lambda_list = logspace(-3,0,10); % List of lambdas to try
    % Emptry matrix to store candidate xi models for this state
    candidate_xi_matrix = zeros(num_functions, length(lambda_list)*size(Theta_i,2)/2+1);
    index = 1; % keeps track of next empty index in candidate_xi_matrix 
    
    for lambda = lambda_list % Test each lambda in list
        tic_lambda=tic();
        for j = 1:1:size(Theta_i,2)/2+1 % Guess a column in Theta (only for numerator)

            Theta_rm = remove_column(Theta_i,j); % Remove column being guessed
            Theta_rm_test = remove_column(Theta_i_test,j);
            
            % sparsifyDynamics hogs the time because of backslash/regression:
            xi = sparsifyDynamics(Theta_rm,Theta_i(:,j),lambda); % Sequential LS
            
            % Calculate model 2-norm error with test data to compare models
            error = norm(Theta_i_test(:,j) - Theta_rm_test*xi, 1)/norm(Theta_i_test(:,j),1);
            
            % Insert -1 into position j for removed column of Theta
            xi_full = [xi(1:j-1); -1; xi(j:end)];
            candidate_xi_matrix(:,index) = xi_full;
            index=index+1;

            hold on;
            % Update best_xi if error is smaller than record
            
            % Cross-validate error bad for online, because of need more data
            % Or maybe you just store previous data, but then slow for
            % changes to dynamics
            % ??? Maybe try change to AIC
            
            if error < min_error
                best_xi = xi_full; % Update xi
                min_error = error; % Update min_error value
                best_lambda = lambda; % Store lambda value used for this column
                best_column = j; % Store column used for this guess
            end
            
            % Calculate Adjusted R^2 value
            % Loosely based on code by; R P (2020). rsquared (https://www.mathworks.com/matlabcentral/fileexchange/60577-rsquared), MATLAB Central File Exchange. Retrieved May 29, 2020.
            
            y_data = Theta_i(:,j); 
            y_estimated = Theta_rm*xi;
            num_param = nnz(xi);
            
            SSres=sum( (y_data - y_estimated).^2 ); % residual sum of squares
            SStot=sum( (y_data - mean(y_data)).^2 ); % total sum of squares
            R2 = 1-SSres/SStot; % standard R squared
            R2adj = 1 - SSres/SStot * (length(y_data)-1)/(length(y_data)-num_param); % adjust for the number of parameters

            if R2adj > max_R2adj
                best_xi_r2 = xi_full; % Update xi
                max_R2adj = R2adj; % Update min_error value
                best_lambda_r2 = lambda; % Store lambda value used for this column
                best_column_r2 = j; % Store column used for this guess
            end
            
        end % End: for each column, j, in Theta
    end % End: for each lambda in lambda_list
    
    % Append model to Xi for this state
    Xi(:,i) = best_xi;
    model_errors(:,i) = min_error;
    model_lambdas(:,i) = best_lambda;
    model_column_guess(:,i) = best_column;
    
    % If adjusted R^2 used as metric
    Xi_r2(:,i) = best_xi;
    model_errors_r2(:,i) = max_R2adj;
    model_lambdas_r2(:,i) = best_lambda_r2;
    model_column_guess_r2(:,i) = best_column_r2;
    
    candidate_xi_cell(i) = {candidate_xi_matrix}; % Add matrix of candidate models for this state to cell
end % End: for each state, i


%% Visualise Xi
x_names = {'x1', 'x2', 'x3', 'sin(x1)'};
% x_names = [x_names, {'1'}, x_names];
vis_Xi = visualize_Xi(x_names, Xi, 2)

model_errors
model_lambdas
model_column_guess % Guessed column of Theta that worked best


%% Run model and compare to actual data

x0 = [-6; 3; 0.8]; % Initial condition for test data
tspan = [0.01:0.01:20];

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t_hat,x_hat] = ode45(@(t,x) SINDy_PI_ODE(t, x, Xi), tspan, x0, options);
[t,x] = ode45(@(t,x) system_ODE(t, x), tspan, x0, options);

figure
plot(t,x); hold on;
plot(t_hat,x_hat,'--', 'LineWidth', 1); hold off;


toc;

warning('on','MATLAB:rankDeficientMatrix'); % Switch on warning for other scripts

function dx = system_ODE(t,x)
    % 2 states
    x1 = x(1,:);
    x2 = x(2,:);
    x3 = x(3,:);
    
    dx = [
            - x1 + 3*x3 - x2.^2;
            (10*x1  - x2 - 2*x1.*x3)./(1 + sin(x1) + x1.^2);
            (x1.*x2 - 3*x3)
    ];

end

function Theta = Theta(X)
%     x1 = X(:,1);
%     x2 = X(:,2);
%     x3 = X(:,3);
%     
%     Theta = [ones(samples,1), x1, x2, x3, x1.*x2, x1.*x3, x1.^2, x2.^2, x3.^2];
    X = [X, sin(X(:,1))];
    Theta = [ones(size(X,1),1), X];
    for i = 1:size(X,2)
        for j = i:size(X,2)
            Theta = [Theta, X(:,i).*X(:,j)];
        end
    end
end

function x_dot = SINDy_PI_ODE(t, x, Xi)
    % ODE created by SINDy_PI model defined by Theta and Xi
    % x is column state vector
    % x_dot is column derivative vector
    n = length(x); % Number of states
    
    X = x'; % Row vector fits SINDY format  
    Theta_X = Theta(X); % Calculate library of functions once
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
    Xi;
    % lambda is our sparsification knob.
    for k=1:5
        small_indexes = (abs(Xi)<lambda);   % find small coefficients
        Xi(small_indexes)=0;                % set small coeffs to 0 (threshold)

        big_indexes = ~small_indexes;
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(big_indexes,:) = Theta(:,big_indexes)\dXdt;
        
        if(Xi == Xi_prev)
            % To save computations:
            % If Xi already converged, then exit loop
            break;
        end
        Xi_prev = Xi; % Save previous Xi

        %         pause
        Xi;
    end
end