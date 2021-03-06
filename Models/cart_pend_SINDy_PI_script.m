%% Performs SINDy-PI on pendulum on a cart system.
% Estimate non-linear model from data
% Requires full state feedback
% Partial state feedback
% Based on code by Eureka Kaiser
% github.com/dynamicslab/SINDy-PI

tic; % Start timer

%% Read Data
load('cartpend_real_Xi'); % Load value for Xi that works
load('cartpend_random_1');

% All simulation data before train/test split
t = out.tout;
Ts = t(2) - t(1);
X = out.x.Data; % All X data (train + test)
X_dot = out.x_dot.Data; % All X_dot data (train + test)
U = out.u.Data; % All U data (train + test)
N = size(X,1); % Number of data samples per state
n = size(X,2); % Number of states

% Parameters
N_train = 2000; % Num of data samples for training, rest for testing
sigma = 0.00000; % Standard deviation of noise

% Train/Test split
X_train = X(1:N_train, :);
X_dot_train = X_dot(1:N_train, :);
U_train = U(1:N_train, :);
t_train = t(1:N_train, :);

N_test = N - N_train + 1; % Num of data samples for testing
X_test = X(N_train:end, :); % One sample overlaps for initial condition
X_dot_test = X_dot(N_train:end, :);
U_test = U(N_train:end, :);
t_test = t(N_train:end, :);

% Add noise to measurements
X_train = X_train + sigma*randn(size(X_train));

%% Total Variation Regularized Differentiation
% Implementation of TVRegDiff from the publication "Sparse identification of nonlinear dynamics for model predictive control in the low-data limit" by E. Kaiser, J. N. Kutz and S. L. Brunton.
denoise = 0; % 1 = Enable denoising
if(denoise)
    alpha = 1*10.^1
    X_dot_clean = zeros(size(X_train)+[1 0]);
    for i = 1:size(X_train,2)
        tic
        X_dot_clean(:,i) = TVRegDiff( X_train(:,i), 10, alpha, [], 'small', 1e6, Ts, 0, 0 ); %.00002
        toc
    end
    % Cutoff index=1, because 'small' setting adds extra entry:
    X_dot_clean = X_dot_clean(2:end,:);

    % Use integral of X_dot_clean for X_clean
    X_clean = zeros(size(X_train));
    for i = 1:size(X_train,2)
        X_clean(:,i) = X_clean(1,i) + cumtrapz(t_train, X_dot_clean(:,i)); % Numeric integration
        X_clean(:,i) = X_clean(:,i) - (mean(X_clean(50:end-50,i)) - mean(X_train(50:end-50,i))); % Adjust mean
    end
    X_clean = X_clean(50:end-51,:);
    X_dot_clean = X_dot_clean(50:end-51,:);  % trim off ends (overly conservative)
    
    tc = t_train(50:end-51);
    for f=1:4
    subplot(1,2,1), plot(t_train,X_dot_train(:,f)); hold on
    plot(tc,X_dot_clean(:,f)); hold off
    title('X dot')
    legend('measured', 'clean')
    subplot(1,2,2), plot(t_train,X_train(:,f)); hold on
    plot(tc,X_clean(:,f)); hold off
    title('X')
    legend('measured', 'clean')
    pause
    end

    % Use denoised data
    X_train = X_clean;
    X_dot_train = X_dot_clean;
    U_train = U_train(50:end-51);
    t_train = t_train(50:end-51);
end

% Plot data
subplot(1,1,1);
figure(1), plot(t_train,X_train);
title("Traing Data");
% drawnow;

%% Find best model for each state

polyorder = 2; % Highest order polynomial term in function library
lambda_list = logspace(-6,-1,5); % List of lambdas to try
% lambda_list = 1e-4;

% Theta to compute function for x_dot
Theta_X = Theta(X_train,U_train,polyorder);
num_functions = size(Theta_X,2)*2; % Number of functions in Theta(x,xdot) library

Xi = zeros(size(Theta_X,2)*2,n); % Stores final Xi, which each column a xi per state
model_errors = []; % error of model for each state
model_lambdas = []; % lambda used for model of each 
model_column_guess = []; % Column used as guess for model of each state

warning('off','MATLAB:rankDeficientMatrix'); % Do not warn about rank deficiency
% k_list = zeros(1,n*length(lambda_list)*(size(Theta_i,2)/2+1));
index = 1;
guess_list = 1:1:size(Theta_X,2);
dont_guess = [24,25]; % Row indices not to guess during regression
guess_list = guess_list(~ismember(guess_list,dont_guess)); % Remove indices of dont_guess

for i = [2, 4] % Only state 1 and 3 % Loop through all states, i
    Theta_i = [Theta_X, diag(X_dot_train(:,i))*Theta_X]; % Theta used for x1 only
    
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
             
%             sparsifyDynamics hogs the time because of backslash/regression:
            [xi,k_conv] = sparsify_dynamics(Theta_rm,Theta_i(:,j),lambda); % Sequential LS
%             k_list(index) = k_conv;
%             xi = Theta_rm\Theta_i(:,j);          
%             full_xi = [xi(1:j-1); -1; xi(j:end)];
%             full_xi = full_xi./sum(full_xi)
%             figure(2), bar3([full_xi, real_Xi(:,i)./sum(real_Xi(:,i))])
%             i
%             lambda
%             j
%             
%             pause
            index = index+1;
            
            % Calculate 2-norm error without test data, penalised for number of parameters
            num_terms = nnz(xi)+1; % Number of non-zero terms in model
            % ??? this error of regression, not of model, can give false
            % low error, like with a trig identity cos^2 + sin^2 = 1
            error = norm(Theta_i(:,j) - Theta_rm*xi, 1)/norm(Theta_i(:,j),1);
            
            metric = error*num_terms^3; % Metric used to compare candidate models
            % Mertric promotes sparsity
            
            % Plot error vs #terms of all candidate models
            
            figure(3), subplot(2,n/2,i/2), semilogy(num_terms,error, color, 'MarkerSize', 8), hold on;
            figure(3), subplot(2,n/2,(i+n)/2), semilogy(num_terms,metric, color, 'MarkerSize', 8), hold on;          

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
    
    figure(3), subplot(2,n/2,i/2), semilogy(nnz(best_xi),best_error1, 'o')
    ylabel('data error');
    grid on;
    
    figure(3), subplot(2,n/2,(i+n)/2), plot(nnz(best_xi),min_metric, 'o')
    ylabel('metric');
    grid on;
    
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
subplot(2,n/2,1), y_scale = ylim;
y_max = max(y_scale);
y_min = min(y_scale);
for i = 2:n
    subplot(2,n/2,i), y_scale = ylim;
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
for i = 1:n
    subplot(2,n/2,i), ylim([y_min y_max]);
end

for i = 1:n
    subplot(2,n/2,i), ylim;
end

%% Compare real Xi to model Xi with bar plots
% load('cartpend_real_Xi')
% index = 2;
% for i = 1:2
%     figure; bar3([Xi(:,index), real_Xi(:,index)]);
%     index = index + 2;
% end

% histogram(k_list) % Display frequency of values of k where Xi converged

%% Visualise Xi
x_names = {'x1', 'x2', 'x3', 'x4', 'u'};
% x_names = [x_names, {'1'}, x_names];

vis_Xi = visualize_Xi(x_names, Xi, polyorder)
vis_real_Xi = visualize_Xi(x_names, real_Xi, polyorder);

lambda_list
model_errors
model_lambdas
model_column_guess % Guessed column of Theta that worked best

disp('Model computation time')
toc; % Display computation time

%% Testing
% Run model on unseen testing data and compare to actual measurements

x0 = X_test(1,:); % Initial condition
x_hat = zeros(N_test,n); % Empty prediction matrix
x_hat(1,:) = x0; 
% Generate data with SINDY-PI model
% Solve for small intervals with constant u
for i=1:N_test-1
    x0 = x_hat(i,:)';
    u = U_test(i,:); %(U_test(i,:) + U_test(i+1,:))/2; % Assume constant u at average of time interval
    [t_1,x_1] = ode45(@(t_1,x_1) SINDy_PI_ODE(t_1,x_1,u,Xi,polyorder), t_test(i:i+1,1), x0);
    x_hat(i+1,:) = x_1(end,:);
end

y_hat = x_hat(:,[1,3]);
y_test = X_test(:,[1,3]);

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat' - y_test'), 2)./N_test % For each measured state

% % Vector of Root Mean Squared Error on testing data
% RMSE = sqrt(sum((y_hat' - y_test').^2, 2)./N_test) % Each row represents RMSE for measured state


%% Plot simulation data vs model data
figure(1), hold on;
plot(t_test,y_test); hold on;
plot(t_test,y_hat,'--', 'LineWidth', 1);
plot(t,U, ':'); 
plot([t(N_train) t(N_train)], ylim, 'k');
hold off;
title('Training and Validation data vs Model');

disp('Total execution time')
toc; % Display total execution time

warning('on','MATLAB:rankDeficientMatrix'); % Switch on warning for other scripts

%% Maximise figures
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
%     X = [X, U, sin(X(:,3)), cos(X(:,3))];
    X = [X, U];
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
        
    % Add extra functions
    Theta_X = [Theta_X, sin(X(:,3)).*cos(X(:,3))];
    Theta_X = [Theta_X, X(:,4).^2.*sin(X(:,3))];
    Theta_X = [Theta_X, cos(X(:,3)).^2];
    Theta_X = [Theta_X, sin(X(:,3))];
    Theta_X = [Theta_X, X(:,4).^2.*sin(X(:,3)).*cos(X(:,3))];
    Theta_X = [Theta_X, X(:,2).*cos(X(:,3))];
    Theta_X = [Theta_X, cos(X(:,3)).*U];
    
    % Remove unneccesary functions
    Theta_X = Theta_X(:,[1:6, 20:end]);
    
    
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

function [Xi,k_conv] = sparsify_dynamics(Theta_X,dXdt,lambda)
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

function vis_Xi = visualize_Xi(x_names, Xi, polyorder)
    % Create cell array to visualise terms in Xi matrix from SINDY
    % x_names   = string names of terms in Theta e.g {'x', 'y', 'z', 'sin(y)', '(1/x)'}
    % Xi        = Xi matrix;
    % polyorder = highest order of polynomial in function

    num_terms = length(x_names);
    vis_Xi = cell(size(Xi)+[1 2]); % Empty visualisation cell matrix
    % Leave space in vis_Xi for 2 columns on left and one on top for headings

    vis_Xi(2:end, 3:end) = num2cell(Xi); % Insert Xi into vis_Xi 
    rows = size(Xi,1); % Number of rows in Xi/ states
    cols = size(Xi,2); % Number of columns in Xi/ terms

    % for i = 1:cols % Add column headers of states
    %     vis_Xi(1,i+2) = strcat(x_names(i),'_dot');
    % end

    vis_Xi(2,2) = {'1'};
    index = 3; % Row index to add next label

    for i=1:num_terms
        vis_Xi(index,2) = x_names(i);
        index = index+1;
    end

    if(polyorder>=2)
        % poly order 2
        for i=1:num_terms
            for j=i:num_terms
                vis_Xi{index,2} = [x_names{i},x_names{j}];
                index = index+1;
            end
        end
    end

    if(polyorder>=3)
        % poly order 3
        for i=1:num_terms
            for j=i:num_terms
                for k=j:num_terms
                    vis_Xi{index,2} = [x_names{i},x_names{j},x_names{k}];
                    index = index+1;
                end
            end
        end
    end
    
    vis_Xi{index,2} = ['(sinx3)(cosx3)'];
    index = index+1;
    vis_Xi{index,2} = ['(x4^2)(sinx3)'];
    index = index+1;
    vis_Xi{index,2} = ['(cosx3)^2'];
    index = index+1;
    vis_Xi{index,2} = ['(sinx3)'];
    index = index+1;
    vis_Xi{index,2} = ['(x4^2)(sinx3)(cosx3)'];
    index = index+1;
    vis_Xi{index,2} = ['(x2)(cosx3)'];
    index = index+1;
    vis_Xi{index,2} = ['(cosx3)(u)'];
    index = index+1;
    
    num_funcs = size(Xi,1);
    
    vis_Xi((floor(end/2)+2):end, 2) = vis_Xi(2:(floor(end/2)+1), 2); % Add same labels from numerator to denomenator

    
    vis_Xi(:,1) = num2cell((0:size(vis_Xi,1)-1)'); % Add row numbering
end
