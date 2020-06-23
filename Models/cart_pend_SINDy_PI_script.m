%% Performs SINDy-PI on pendulum on a cart system.
% Based on code by Eureka Kaiser
% github.com/dynamicslab/SINDy-PI

%% First try with simple system called "system_ODE"
% close all;
clear all;

tic; % Start timer

%% Read Data

% Load: x0 and out from simulation

% load('rational_toy_with_input_data_1.mat') % Polyorder = 2
% load('rational_toy_poly3_1.mat') % Polyorder = 3
load('cartpend_data_3')
load('cartpend_real_Xi') % Load value for Xi that works

t = out.tout;
Ts = t(2) - t(1);
X = out.x.Data;
X_dot = out.x_dot.Data; % ??? Change to calculate dx with total variation derivative
U = out.u.Data;
num_samples = size(X,1); % Number of data samples per state
n = size(X,2); % Number of states

% Add noise to measurements
sigma   = 0.001; % Magnitude of noise
X       = X + sigma*randn(size(X));

%% Total Variation Regularized Differentiation
% Implementation of TVRegDiff from the publication "Sparse identification of nonlinear dynamics for model predictive control in the low-data limit" by E. Kaiser, J. N. Kutz and S. L. Brunton.
denoise = 1;
if(denoise)
    alpha = 7*10.^1
    X_dot_clean = zeros(size(X)+[1 0]);
    for i = 1:size(X,2)
        tic
        X_dot_clean(:,i) = TVRegDiff( X(:,i), 10, alpha, [], 'small', 1e6, Ts, 0, 0 ); %.00002
        toc
    end
    % Because 'small' adds extra entry:
    X_dot_clean = X_dot_clean(2:end,:);

    % Use integral of X_dot_clean for X_clean
    X_clean = zeros(size(X));
    for i = 1:size(X,2)
        X_clean(:,i) = X_clean(1,i) + cumtrapz(t, X_dot_clean(:,i)); % Numeric integration
        X_clean(:,i) = X_clean(:,i) - (mean(X_clean(50:end-50,i)) - mean(X(50:end-50,i))); % Adjust mean
    end
    X_clean = X_clean(50:end-51,:);
    X_dot_clean = X_dot_clean(50:end-51,:);  % trim off ends (overly conservative)
    % 
    % tc = t(50:end-51);
    % f=1;
    % for f=1:4
    % subplot(1,2,1), plot(t,X_dot(:,f)); hold on
    % plot(tc,X_dot_clean(:,f)); hold off
    % title('X dot')
    % legend('measured', 'clean')
    % subplot(1,2,2), plot(t,X(:,f)); hold on
    % plot(tc,X_clean(:,f)); hold off
    % title('X')
    % legend('measured', 'clean')
    % pause
    % end

    % Use denoised data
    X = X_clean;
    X_dot = X_dot_clean;
    U = U(50:end-51);
    t = t(50:end-51);
end


% Choose window size for training data
window = 5000;
X = X(1:window,:);
X_dot = X_dot(1:window,:);
U = U(1:window,:);
t = t(1:window,:);

% Plot data
figure(1), plot(t,X);
title("Training data");
drawnow;

%% Find best model for each state

polyorder = 2; % Highest order polynomial term in function library
lambda_list = logspace(-3,-1,4); % List of lambdas to try
% lambda_list = 1e-4;

% Theta to compute function for x_dot
Theta_X = Theta(X,U,polyorder);
num_functions = size(Theta_X,2)*2; % Number of functions in Theta(x,xdot) library

Xi = zeros(size(Theta_X,2)*2,n); % Stores final Xi, which each column a xi per state
model_errors = []; % error of model for each state
model_lambdas = []; % lambda used for model of each 
model_column_guess = []; % Column used as guess for model of each state

warning('off','MATLAB:rankDeficientMatrix'); % Do not warn about rank deficiency
% k_list = zeros(1,n*length(lambda_list)*(size(Theta_i,2)/2+1));
index = 1;
guess_list = 1:1:size(Theta_X,2);
guess_list = guess_list(guess_list~=1);
guess_list = guess_list(guess_list~=34); % Remove sin^2 because trig identity give false low error
guess_list = guess_list(guess_list~=35);
guess_list = guess_list(guess_list~=36); % Remove cos^2
guess_list = [3,6, 7];

for i = [2, 4] % Only state 1 and 3 % Loop through all states, i
    Theta_i = [Theta_X, diag(X_dot(:,i))*Theta_X]; % Theta used for x1 only
    
    min_metric = Inf; % Store minimum model 2-norm error
    best_lambda = Inf; % Store best lambda
    best_column = Inf; % Store best guess of column
    best_xi = []; % Store best xi for model
   
    color_list = {'b.', 'g.', 'k.', 'm.', 'r.', 'y.'...
        'bx', 'gx', 'kx', 'mx', 'rx', 'yx'...
        'b+', 'g+', 'k+', 'm+', 'r+', 'y+'}; % Range of colors to use
    
    color_index = 1; % Index of next color to use
    for lambda = lambda_list % Test each lambda in list
        color = color_list{color_index}; % New color for each lambda
        color_index = color_index + 1; % Next color index
        tic_lambda = tic();
        for j = guess_list % Guess a column in Theta (only for numerator)

            Theta_rm = remove_column(Theta_i,j); % Remove column being guessed
             
            % sparsifyDynamics hogs the time because of backslash/regression:
%             [xi,k_conv] = sparsifyDynamics(Theta_rm,Theta_i(:,j),lambda); % Sequential LS
%             k_list(index) = k_conv;
            xi = Theta_rm\Theta_i(:,j);          
            full_xi = [xi(1:j-1); -1; xi(j:end)];
            full_xi = full_xi./sum(full_xi)
            figure(2), bar3([full_xi, real_Xi(:,i)./sum(real_Xi(:,i))])
            i
            lambda
            j
            
            pause
            index = index+1;
            
            % Calculate 2-norm error without test data, penalised for number of parameters
            num_terms = nnz(xi)+1; % Number of non-zero terms in model
            % ??? this error of regression, not of model, can give false
            % low error, like with a trig identity cos^2 + sin^2 = 1
            error = norm(Theta_i(:,j) - Theta_rm*xi, 1)/norm(Theta_i(:,j),1);
            
            metric = error*num_terms^3; % Metric used to compare candidate models
            % Mertric promotes sparsity
            
            % Plot error vs #terms of all candidate models
            
            figure(3), subplot(2,n,i), semilogy(num_terms,error, color, 'MarkerSize', 8), hold on;
            figure(3), subplot(2,n,i+n), semilogy(num_terms,metric, color, 'MarkerSize', 8), hold on;          

            % ??? Maybe try change to AIC
            
            % Update best_xi if metric is smaller than smallest metric yet
            if metric < min_metric
                best_xi = [xi(1:j-1); -1; xi(j:end)];  % Insert -1 into position j for removed column of Theta
                min_metric = metric; % Update min_error value
                best_error1 = error; % Save 2-norm error of original data
                best_lambda = lambda; % Store lambda value used for this column
                best_column = j; % Store column used for this guess
            end
            
        end % End: for each column, j, in Theta
        
    end % End: for each lambda in lambda_list
    
    figure(3), subplot(2,n/2,i/2), plot(nnz(best_xi),best_error1, 'o')
    ylabel('data error');
    grid on;
    hold off;
    
    figure(3), subplot(2,n/2,(i+n)/2), plot(nnz(best_xi),min_metric, 'o')
    ylabel('metric');
    grid on;
    hold off;
    
    % Append model to Xi for this state
    Xi(:,i) = best_xi;
    model_errors(:,i) = min_metric;
    model_lambdas(:,i) = best_lambda;
    model_column_guess(:,i) = best_column;
    
end % End: for each state, i

%%
% Manually add states 1 and 3 
Xi(3,1) = -1; % x1_dot = x2
Xi(5,3) = -1; % x3_dot = x4

Xi(end/2+1,1) = 1; % Add denominator = 1
Xi(end/2+1,3) = 1;

% If SINDy-PI added no terms to denominator, add them
for i = 1:n
    if nnz(Xi((end/2+1):end,i)) == 0 % If no terms in denominator
        Xi((end/2+1),i) = 1; % Add denominator = 1
    end
end

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
    subplot(2,n,i), ylim;
end

%% Compare real Xi to model Xi with bar plots
load('cartpend_real_Xi')
figure;
index = 2;
for i = 1:2
    subplot(1,2,i), bar3([Xi(:,index), real_Xi(:,index)]);
    index = index + 2;
end

% histogram(k_list) % Display frequency of values of k where Xi converged

%% Visualise Xi
x_names = {'x1', 'x2', 'x3', 'x4', 'u', 'sin(x3)', 'cos(x3)'};
% x_names = [x_names, {'1'}, x_names];

vis_Xi = visualize_Xi(x_names, Xi, polyorder)

lambda_list
model_errors
model_lambdas
model_column_guess % Guessed column of Theta that worked best

disp('Model computation time')
toc; % Display computation time

%% Validatation
% Run model on new data and compare to actual measurements

% Load x0, out from a simulation
% load('rational_toy_with_input_data_2.mat') % Polyorder = 2
% load('rational_toy_poly3_2.mat') % Polyorder = 3
load('cartpend_data_4')

t_valid = out.tout;
X_valid = out.x.Data;
U_valid = out.u.Data;

% Only use portion of data
t_valid = t_valid(1:window,:);
X_valid = X_valid(1:window,:);
U_valid = U_valid(1:window,:);

% Generate data with SINDY-PI model
x0 = X_valid(1,:);
x_hat = zeros(window,n);
t_hat = zeros(window,1);
x_hat(1,:) = x0;
t_hat(1,:) = 0;

% Solve for small intervals with constant u
for i=2:size(U_valid,1)
    u = (U_valid(i-1,:)+U_valid(i,:))/2; % Assume constant u at average of time interval
    [t_1,x_1] = ode45(@(t_1,x_1) SINDy_PI_ODE(t_1,x_1,u,Xi,polyorder), t_valid(i-1:i,1), x0);
    x_hat(i,:) = x_1(end,:);
    t_hat(i,:) = t_1(end,:);
    x0 = x_hat(i,:)';
end

% PLot simulation data vs model data
figure
plot(t_valid,X_valid); hold on;
plot(t_hat,x_hat,'--', 'LineWidth', 1); hold off;
title('Validation data vs Model');

disp('Total execution time')
toc; % Display total execution time

warning('on','MATLAB:rankDeficientMatrix'); % Switch on warning for other scripts

% Maximise figures
for i=1:4
    figure(i), set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

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

function Theta_X = Theta(X, U, polyorder)
    %  Poly order =  Highest order of polynomial term in library
    %   Theta = [ones(samples,1), x1, x2, x3, x1.*x2, x1.*x3, x1.^2, x2.^2, x3.^2];
    X = [X, U, sin(X(:,3)), cos(X(:,3))];
    
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
    
    Theta_X = [Theta_X, X(:,4).^2.*sin(X(:,3)), X(:,4).^2.*sin(X(:,3)).*cos(X(:,3))];
    
end

function x_dot = SINDy_PI_ODE(t, x, u, Xi, polyorder)
    % ODE created by SINDy_PI model defined by Theta and Xi
    % x is column state vector
    % x_dot is column derivative vector
    n = length(x); % Number of states
    
    X = x'; % Row vector fits SINDY format  
    Theta_X = Theta(X,u,polyorder); % Calculate library of functions once
    X_dot = zeros(size(X)); % Allocate space once
    
    for i = 1:n % For each state
        xi = Xi(:,i); 
        xi_num = xi(1:(end/2)); % Top half of xi for model numerator
        xi_den = xi((end/2+1):end); % Bottom half of xi for model denomenator
  
        X_dot(:,i) = -(Theta_X*xi_num)./(Theta_X*xi_den);
    end
    
    x_dot = X_dot'; % Back to column vector
end

function Theta_rm = remove_column(Theta_X,column)
    % Remove the column of Theta given by index: column
    num_columns = size(Theta_X,2);
    Theta_rm = Theta_X(:,(1:num_columns)~=column);
end

function [Xi,k_conv] = sparsifyDynamics(Theta_X,dXdt,lambda)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz

    % compute Sparse regression: sequential least squares
    Xi = Theta_X\dXdt;  % initial guess: Least-squares
    Xi_prev = Xi;
    
    % lambda is our sparsification knob.
    k_conv = 20; % k at which Xi converges
    for k=1:k_conv
        % lambda is a fraction of the largest value in the model which
        % is considered significant
        % Change from original lambda
        threshold = lambda*max(Xi);
        small_indexes = (abs(Xi)<threshold);   % find small coefficients
        Xi(small_indexes)=0;                % set small coeffs to 0 (threshold)

        big_indexes = ~small_indexes;
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(big_indexes,:) = Theta_X(:,big_indexes)\dXdt;
        
        if(Xi == Xi_prev)
            % To save computations, almost by half:
            % If Xi already converged, then exit loop
            k_conv = k;
            break;
        end
        Xi_prev = Xi; % Save previous Xi
    end
end


