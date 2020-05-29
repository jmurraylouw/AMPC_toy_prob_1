%% Performs SINDy-PI on pendulum on a cart system.
% Based on code by Eureka Kaiser
% github.com/dynamicslab/SINDy-PI

%% First try with simple system called "system_ODE"
tic;
%% Generate Data
x0 = [0.1; -0.2; 0.3];  % Initial condition
x0_test = [1; 0.1; -0.5]; % Initial condition for test data

n = length(x0);  % Number of states
tspan = [0.01:0.01:10];

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x] = ode45(@(t,x) system_ODE(t,x), tspan, x0, options);
[t_test,x_test] = ode45(@(t,x) system_ODE(t,x), tspan, x0_test, options);

subplot(1,2,1), plot(t,x);
subplot(1,2,2), plot(t_test,x_test);

% Calculate derivatives
dx = system_ODE(t,x');
dx_test = system_ODE(t_test,x_test');

% State and derivative matrixes
X = x;
X_dot = dx';

X_test = x_test;
X_dot_test = dx_test';

% Theta to compute function for x_dot
Theta_X = Theta(X);
Theta_X_test = Theta(X_test);

%% Find best model for each state
Xi = []; % Stores final Xi, which each column a xi per state
model_errors = []; % error of model for each state
model_lambdas = []; % lambda used for model of each 
model_column_guess = []; % Column used as guess for model of each state

for i = 1:n
    Theta_i = [Theta_X, diag(X_dot(:,i))*Theta_X]; % Theta used for x1 only
    Theta_i_test = [Theta_X_test, diag(X_dot_test(:,i))*Theta_X_test]; % test Theta used for x1 only
    
    min_error = Inf; % Store minimum model 2-norm error
    best_lambda = Inf; % Store best lambda
    best_column = Inf; % Store best guess of column
    best_xi = []; % Store best xi for model

    for lambda = 1e-7%:1e-7:1e-5 % Test each lambda in list
        tic_lambda=tic();
        for j = 1:1:size(Theta_i,2) % Guess every column in Theta
            Theta_rm = remove_column(Theta_i,j); % Remove column being guessed
            Theta_rm_test = remove_column(Theta_i_test,j);

    %         xi = lasso(Theta_rm,Theta_1(:,j),'Lambda',lambda);
            xi = sparsifyDynamics(Theta_rm,Theta_i(:,j),lambda,1); % Sequential LS

            % Calculate model 2-norm error with test data to compare models
            error = norm(Theta_i_test(:,j) - Theta_rm_test*xi, 1)/norm(Theta_i_test(:,j),1);

            % Update best_xi if error is smallest yet
            if error < min_error
                % Insert -1 into position j for removed column of Theta
                best_xi = [xi(1:j-1); -1; xi(j:end)];
                min_error = error; % Update min_error value
                best_lambda = lambda; % Store lambda value used for this column
                best_column = j; % Store column used for this guess
            end
        end % End: for each column, j, in Theta
    end % End: for each lambda in lambda_list
    
    % Append model to Xi for this state
    Xi(:,i) = best_xi;
    model_errors(:,i) = min_error;
    model_lambdas(:,i) = best_lambda;
    model_column_guess(:,i) = best_column;

end % End: for each state, i

%% Visualise Xi
x_names = {'x1', 'x2', 'x3', 'x1.*x2', 'x1.*x3', 'x1.^2', 'x2.^2', 'x3.^2'};
x_names = [x_names, {'1'}, x_names];
vis_Xi = visualize_Xi(x_names, Xi, 1)

model_errors
model_lambdas
x_names(model_column_guess-1) % Guessed column of Theta that worked best


%% Run model and compare to actual data

x0 = [-6; 3; 0.8]; % Initial condition for test data
tspan = [0.01:0.01:20];

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x] = ode45(@(t,x) system_ODE(t, x), tspan, x0, options);
[t_hat,x_hat] = ode45(@(t,x) SINDy_PI_ODE(t, x, Xi), tspan, x0);

figure
plot(t,x); hold on;
plot(t_hat,x_hat,'--'); hold off;


toc;

function dx = system_ODE(t,x)
    % 2 states
    x1 = x(1,:);
    x2 = x(2,:);
    x3 = x(3,:);
    
    dx = [
            - x1 + 3*x3 - x2.^2;
            (10*x1  - x2 - 2*x1.*x3)./(1+x1.^2);
            x1.*x2 - 3*x3
    ];

end

function Theta = Theta(X)
    x1 = X(:,1);
    x2 = X(:,2);
    x3 = X(:,3);
    samples = length(x1); % Number of samples per state (number of rows)
    
    Theta = [ones(samples,1), x1, x2, x3, x1.*x2, x1.*x3, x1.^2, x2.^2, x3.^2];
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
