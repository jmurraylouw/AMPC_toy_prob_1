%% DMD - Moving-Window of cart pendulum
% Estimate System matrixes with moving window of data in real time
% Full state feedback

u_data  = out.u.Data';
x_data  = out.x.Data';
t       = out.tout';

%%
N       = max(size(x_data));
Ts      = 0.05;     % Sample time of data
w       = 100;    % Size of window in timesteps

%%
% Initialise
A_output = zeros(2,2);
B_output = zeros(2,1);

X_dwork = zeros(2,w);
U_dwork = zeros(1,w);

% plot(t, x_data, 'k'); hold on; % Plot measured x vs t

a11=zeros(1,N-1);
MSE = zeros(1,N-1);

incr = 1;%/Ts; % Increment size for data plots
for k = 1:incr:N-1
    % Inport Dwork memory
    X = X_dwork;
    U = U_dwork;
    
    % Inputs
    x       = x_data(k);
    u       = u_data(k);
    
    % Add input data to X2
    X2 = [X(:, 2:end), [x_dot; x]];
   
    % Based on DMD control example video by Steve Brunton
    XU = [X; U];
    AB = X2*pinv(XU);
    A  = AB(:,1:2);
    B  = AB(:,end);
    
    % Check for change in model
    w_c = 20; % window of data points to check for change in system model
    
    % Take only last w_c entries of X and X2
    X2_measured = X2(:, (end-w_c+1):end);
    X_measured  = X(:, (end-w_c+1):end);
    
    % Calculate X2 according to A
    X2_calc = A*X_measured;
    MSE(k) = mean((X2_measured - X2_calc).^2, 'all');
    
    a11(k) = A(1,1);
%     plot_model(A,B,U_data,t);
%     k
%     A
%     B
%     pause
    
    % Output
    A_output = A;
    B_output = B;
    
    % Update Dwork memory
    X_dwork = X2;
    U_dwork = [U(:, 2:end), u];
    
end

C = [0 1];
D = 0;

t = 1:1:N-1;
plot(t,a11);
hold on;
plot(t,MSE)

% plot_model(A,B,f_data,t,2)
% hold on;
% [A,B,C,D] = ssdata(mpc1.Model.Plant)
% plot_model(A,B,U_data,t,2)


A = [0.969169519504925  -0.246386829627017;  0.049277365925403  0.993808202467626];
B = [0.049277365925403; 0.001238359506475];

%%

% Extract nx, ny, nu
size_B = size(B);
nx = size_B(1);
nu = size_B(2);
size_C = size(C);
ny = size_C(1);

% Create nominal point at all 0, because linear model.
X = zeros(nx,1);
Y = zeros(ny,1);
U = zeros(nu,1);
DX = [0; 0];

%% Local functions
function plot_model(A,B,U_data,t,r)
    n = max(size(t));   % Number of time steps 
    X_hat = zeros(2,n); % Estimated X from model
    X_hat(:,1) = [0; 0]; % Initial conditions
    for index = 1:1:n-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
    x_hat = X_hat(r,:); % First row is x
    plot(t, x_hat);     % Estimated x
end







