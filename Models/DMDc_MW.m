%% Least Squares - Moving-Window
% Estimate System matrixes with moving window of data in real time
% Similar to DMD method
% Full state feedback
% To be upgraded to sfunc to use in simulink

f_data      = out.f.Data';
x_data      = out.x.Data';
x_dot_data  = out.x_dot.Data';
t       = out.tout';

%%
n       = max(size(x_data));
Ts      = t(2)-t(1); % Sample time of data
Ts_c    = 0.1;
T_window = 1;               % Size of window in time
w       = (T_window/Ts);    % Size of window in timesteps

% Initialise
A_output = zeros(2,2);
B_output = zeros(2,1);

X_dwork = zeros(2,w);
U_dwork = zeros(1,w);

plot(t, x_data, 'k'); hold on; % Plot measured x vs t

incr = 1;%/Ts; % Increment size for data plots
for k = 1:incr:n-1
    % Inport Dwork memory
    X = X_dwork;
    U = U_dwork;
    
    % Inputs
    x       = x_data(k);
    x_dot   = x_dot_data(k);
    f       = f_data(k);
    
    % Add input data to X2
    X2 = [X(:, 2:end), [x_dot; x]];
   
    % Based on DMD control example video by Steve Brunton
    XU = [X; U];
    AB = X2*pinv(XU);
    A  = AB(:,1:2);
    B  = AB(:,end);
    
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
    U_dwork = [U(:, 2:end), f];
    
end

C = [0 1];
D = 0;

plot_model(A,B,f_data,t,2)
hold on;
% [A,B,C,D] = ssdata(mpc1.Model.Plant)
% plot_model(A,B,U_data,t,2)

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







