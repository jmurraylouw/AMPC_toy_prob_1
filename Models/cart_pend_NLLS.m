% Estimates parameters for cartpend using nlgreyest method from MATLAB

% Setup
rng(0);
random_generator = rng; % repeatable random numbers

% Read simulation data
load('cartpend_random_1');
x_data = out.x.Data'; % each row is timeseries of a state
y_data = x_data([1,2,3,4],:);
u_data = out.u.Data';
t = out.x.Time';
Ts = t(2)-t(1);
N = length(t); % Total number of samples in data

% Testing data - Last 50 s is for testing and one sample overlaps training
N_test = 5000; % Num of data samples for testing
% One sample of testing data overlaps for initial condition
x_test = x_data(:,end-N_test+1:end)'; % each row is timeseries of a state
y_test = y_data(:,end-N_test+1:end)'; 
u_test = u_data(:,end-N_test+1:end)';
t_test = t(:,end-N_test+1:end)';

% Training data - Last sample of training is first sample of testing
N_train = 2000; % Num of data samples for training, rest for testing
x_train = x_data(:,end-N_test-N_train+2:end-N_test+1)';
y_train = y_data(:,end-N_test-N_train+2:end-N_test+1)';
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1)';
t_train = t(:,end-N_test-N_train+2:end-N_test+1)';

% Add noise
sigma = 0.01;
y_train = y_train + sigma*randn(size(y_train));

% Actaul model parameter values
M = 4;
m = 2;
L = 1;
d = 5;
g = 9.81;

P_actual_2 = [  L*m;
                1;
                -d;
                g*m;
                m;
                -(M + m)]./((M + m));
                
P_actual_4 = [  -L*m;
                -1;
                d;
                -(g*m + M*g);
                L*m ;
                -(L*(M + m))]./(L*(M + m));

P_actual = zeros(6,4);
P_actual(:,2) = P_actual_2;
P_actual(:,4) = P_actual_4;

%% Total Variation Regularized Differentiation
% Implementation of TVRegDiff from the publication "Sparse identification of nonlinear dynamics for model predictive control in the low-data limit" by E. Kaiser, J. N. Kutz and S. L. Brunton.

denoise = 1; % 1 = Enable denoising
if(denoise)
    x_train = y_train; % Noisy data
    
    alpha = 1*10.^-1; % works for sigma = 0.01 and 0.1
    X_dot_clean = zeros(size(x_train)+[1 0]);
    for i = 1:size(x_train,2)
        tic
        X_dot_clean(:,i) = TVRegDiff( x_train(:,i), 10, alpha, [], 'small', 1e6, Ts, 0, 0 ); %.00002
        toc
    end
    % Cutoff index=1, because 'small' setting adds extra entry:
    X_dot_clean = X_dot_clean(2:end,:);

    % Use integral of X_dot_clean for X_clean
    X_clean = zeros(size(x_train));
    for i = 1:size(x_train,2)
        X_clean(:,i) = X_clean(1,i) + cumtrapz(t_train, X_dot_clean(:,i)); % Numeric integration
        X_clean(:,i) = X_clean(:,i) - (mean(X_clean(50:end-50,i)) - mean(x_train(50:end-50,i))); % Adjust mean
    end
    X_clean = X_clean(50:end-51,:);
    X_dot_clean = X_dot_clean(50:end-51,:);  % trim off ends (overly conservative)
    
    tc = t_train(50:end-51);
%     for f=1:4
%         subplot(1,2,1), plot(t_train,X_dot_train(:,f)); hold on
%         plot(tc,X_dot_clean(:,f)); hold off
%         title('X dot')
%         legend('measured', 'clean')
%         subplot(1,2,2), plot(t_train,X_train(:,f)); hold on
%         plot(tc,X_clean(:,f)); hold off
%         title('X')
%         legend('measured', 'clean')
%         pause
%     end

    % Use denoised data
    x_train = X_clean;
    x_dot_train = X_dot_clean;
    y_train = y_train(50:end-51,:);
    u_train = u_train(50:end-51,:);
    t_train = t_train(50:end-51,:);
else
    % Calculate x_dot with actual model (Can also export from simulation)
    x_train = y_train; % Noisy data
    x_dot_train = estimated_ODE(t_train',x_train',u_train',P_actual)';
end

% Plot data
figure(1);
plot(t_train, y_train);
hold on;
plot(t_train, x_train);
hold off;
title("Training Data");

%% Non-Linear Least Squares for rational function
% Estimate coeffs for terms in ode

% Based on paper: 
% Recursive Parameter Estimation for Nonlinear Rational Models.
% Online URL for this paper:
% http://eprints.whiterose.ac.uk/78722/
% Q.M. Zhu
Theta = zeros(5,1); % Parameter vector
for state = [2,4]
    switch state
        case 2           
            % Numerator functions
            p_n = [   sin(x_train(:,3)).*x_train(:,4).^2,...    % Th_n(1) = matching param
                      u_train(:,1),...                          % Th_n(2)
                      x_train(:,2),...                          % Th_n(3)
                      cos(x_train(:,3)).*sin(x_train(:,3)) ];   % Th_n(4)

            % Denomenator functions
            p_d = [   cos(x_train(:,3)).^2,...    % Th_d(1) = matching para
                      ones(size(x_train(:,3))) ]; % Th_d(2)  
        case 4
            % Numerator functions
            p_n = [   cos(x_train(:,3)).*sin(x_train(:,3)).*x_train(:,4).^2,...    % Th_n(1) = matching param
                      cos(x_train(:,3)).*u_train(:,1),...                          % Th_n(2)
                      cos(x_train(:,3)).*x_train(:,2),...                          % Th_n(3)
                      sin(x_train(:,3)) ];   % Th_n(4)

            % Denomenator functions
            p_d = [   ones(size(x_train(:,3))),...    % Th_d(1) = matching para
                      cos(x_train(:,3)).^2  ]; % Th_d(2)
    end
    
    % Phi
    Phi_n = p_n;
    Phi_d = p_d.*x_dot_train(:,state);
%     Phi = [Phi_n, Phi_d];
    Phi_d1 = Phi_d(:,1);        % 1st column of Phi_d
    Phi_dr = Phi_d(:,2:end);    % Rest of Phi_d
    
    Theta(:,state) = [p_n, Phi_dr]\Phi_d1;
    
    % Psi and psi
    num = size(p_n,2);
    den = size(p_d,2);
    
   

    % Least squares 
    
%     Y = x_dot_train(:,state);
%     Theta(:,state) = inv(Phi'*Phi)*(Phi'*Y);
    
end
Theta

% Actual values
Th_n2 = [  L*m;
           1;
           -d;
           g*m ];
       
Th_d2 = [ -m;
          (M + m) ];
      
Theta2 = [Th_n2; Th_d2];
Theta2

%% Testing
% Run model on unseen testing data and compare to actual measurements

disp('Testing...')
tic;
x0 = x_test(1,:); % Initial condition
nx = length(x0);
x_hat = zeros(N_test,nx); % Empty prediction matrix
x_hat(1,:) = x0; 

% Generate data with estimated model
% Solve for small intervals with constant u
for i=1:N_test-1
    x0 = x_hat(i,:)';
    u = u_test(i,:); %(U_test(i,:) + U_test(i+1,:))/2; % Assume constant u at average of time interval
    [t_1,x_1] = ode15s(@(t_1,x_1) estimated_ODE(t_1,x_1,u,Theta), t_test(i:i+1,1), x0);
%     [t_1,x_1] = ode45(@(t_1,x_1) ode(t_1,x_1,u), t_test(i:i+1,1), x0);
    x_hat(i+1,:) = x_1(end,:);
end

toc;
%%
y_hat = x_hat(:,[1,3]);
y_test = x_test(:,[1,3]);

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat' - y_test'), 2)./N_test % For each measured state
sigma
%% Plot simulation data vs model data
figure(2), hold on;
plot(t_test,y_test); hold on;
plot(t_test,y_hat,'--', 'LineWidth', 1);
% plot(t_test,u_test, ':'); 
plot([t(N-N_test-N_train) t(N-N_test-N_train)], ylim, 'r');
plot([t(N-N_test) t(N-N_test)], ylim, 'k');
hold off;
title('Training and Validation data vs NLLS Model');


%% Functions
function dx = estimated_ODE(t,x,u,Theta)
    % ODE for cartpend model with estimated parameters
    N = size(x,2);
    dx = zeros(4,N);
    
    dx(1,:) = x(2,:);
    dx(2,:) = ( Theta(1,2).*sin(x(3,:)).*x(4,:).^2              + Theta(2,2).*u(1,:)              + Theta(3,2).*x(2,:)              + Theta(4,2).*cos(x(3,:)).*sin(x(3,:)) )./(  Theta(5,2).*cos(x(3,:)).^2 - 1);
    dx(3,:) = x(4,:);
    dx(4,:) = ( Theta(1,4).*cos(x(3,:)).*sin(x(3,:)).*x(4,:).^2 + Theta(2,4).*u(1,:).*cos(x(3,:)) + Theta(3,4).*x(2,:).*cos(x(3,:)) + Theta(4,4).*sin(x(3,:))              )./(  Theta(5,4).*cos(x(3,:)).^2 - 1);

end

function dx = actual_ode(t,x,u)
    %CARTPEND Models a continuous system of a pendulem on a cart.
    %   based on Steve Brunton code. See youtube.com/watch?v=qjhAAQexzLg&list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m&index=12
    %   x  = state vector [x; x_dot; theta; theta_dot]
    %   dx = derivative of state vector
    %   u  = input vector [f]
    %   m  = mass of pendulem end
    %   M  = mass of cart
    %   L  = length of pendulem rod
    %   g  = acceleration due to gravity
    %   d  = damping coef of friction on cart
    %,m,M,L,g,d,u
    
    M = 4;
    m = 2;
    L = 1;
    d = 5;
    g = 9.81;

    
    % Derivatives
    dx = zeros(4,1);

    Sx = sin(x(3));
    Cx = cos(x(3));
    
    % Equations from derive_cartpend.m
    dx(1,1) = x(2);
    dx(2,1) = (L*m*sin(x(3))*x(4)^2 + u - d*x(2) + g*m*cos(x(3))*sin(x(3)))/(- m*cos(x(3))^2 + M + m);
    dx(3,1) = x(4);
    dx(4,1) = -(L*m*cos(x(3))*sin(x(3))*x(4)^2 + u*cos(x(3)) - d*x(2)*cos(x(3)) + g*m*sin(x(3)) + M*g*sin(x(3)))/(L*(- m*cos(x(3))^2 + M + m));
end







