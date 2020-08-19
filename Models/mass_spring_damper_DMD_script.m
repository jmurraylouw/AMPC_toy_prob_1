%% DMD - Moving-Window of mass spring damper
% Estimate System matrixes with moving window of data in real time
% Partial state feedback

%% Variables for simulation
x0 = [0; 0];
m = 1;
b = 0.1;
k = 5;

%% Read data
% u_data  = out.u.Data';
% x_data  = out.x.Data';
% y_data  = [1 0]*out.x.Data';
% % y_data  = out.y.Data';
% t       = out.tout';

load('C:\Users\Murray\OneDrive - Stellenbosch University\Masters\AMPC_toy_prob_1\Data\mass_spring_damper_Square_wave_input_Measure_position.mat')

nx = size(x_data)*[1; 0];
ny = size(y_data)*[1; 0];
nu = size(u_data)*[1; 0];

N       = max(size(x_data));
Ts      = t(2)-t(1);     % Sample time of data
w       = 100;    % Size of window in timesteps


% %% Batch DMD - Full state feedback, no truncation
% X = x_data(:, 1:end-1);
% X2 = x_data(:, 2:end);
% Upsilon = u_data(:, 1:end-1); % Upsilon
% % X2 = A*X + B*U
% % X2 = [A, B]*[X; U]
% 
% G = X2/[X; Upsilon];
% A = G(:, 1:nx);
% B = G(:, nx+1:end);
% 
% x_hat_data = plot_model(A,B,u_data,t,x0);
% hold on;
% plot(t, x_data, '--')
% hold on;
% 
% MSE_dmd = mean(((x_data-x_hat_data).^2)')'

%% Batch DMDc - Partial state feedback
% Augment y with time delay coordinates of y

samples = 1000;
delays = 3; % Number of delay cordinates, including y_data(1:samples+1)
tau = 1; % Sample number shift of delay cordinates.
% Noticed error increased for increasing tau
% i.e Y = y(k)      y(k+1)      y(k+2)...
%         y(k+tau)  y(k+tau+1)  y(k+tau+2)...

assert(samples+delays*tau <= N); %% otherwise index out of bounds 

% Step 1: Collect and construct the snapshot matrices:
H = []; % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
for i = 1:tau:delays*tau
    H = [y_data(:, i:samples+i); H];  
end

n = size(H,1); % Length of augmented state vector

X2 = H(:, 2:end); % Y advanced 1 step into future
X = H(:, 1:end-1); % Cut off last sample

Upsilon = u_data(:, (delays*tau):(samples + delays*tau - 1)); % Upsilon

Omega = [X; Upsilon]; % Omega is concatination of Y and Upsilon

% Step 2: Compute and truncate the SVD of the input space Omega
[U,S,V] = svd(Omega, 'econ');

% plot(diag(S), 'o')
p = size(V)*[0; 1]; % Reduced rank of Omega svd
U_tilde = U(:, 1:p); % Truncate SVD matrixes of Omega
S_tilde = S(1:p, 1:p);
V_tilde = V(:, 1:p);
U1_tilde = U_tilde(1:n, :);
U2_tilde = U_tilde(n+1:end, :);

% Step 3: Compute the SVD of the output space X'
[U,S,V] = svd(X2, 'econ');

% plot(diag(S), 'o')
r = p-1; % Reduced rank of X2 svd, r < p
U_hat = U(:, 1:r); % Truncate SVD matrixes of X2
S_hat = S(1:r, 1:r);
V_hat = V(:, 1:r);

% Step 4: Compute the approximation of the operators G = [A B]
A_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U1_tilde'*U_hat;
B_tilde = U_hat'*X2*V_tilde*inv(S_tilde)*U2_tilde';

% x_tilde(k+1) = A_tilde*x_tilde(k) + B_tilde*u(k)
% x = U_hat*x_tilde, Transform to original coordinates
% x_tilde = U_hat'*x, Transform to reduced order coordinates
% Here x is augmented state

A = U_hat*A_tilde*U_hat'
B = U_hat*B_tilde

% Manually Add row for x_dot, with Centred Euler differentiation
x_dot_row = 1/(2*Ts)*[1 0 0]*(A^2-eye(3));
A = [x_dot_row; A];
A = [zeros(4,1), A]; % Add column of zeros, nothing depends on x_dot
B = [0; B]; % Add zero. Force does not affect x_dot in model

% Ignore eigenmodes Step 5 and 6

% Run model
x_hat_0 = [];
for i = 1:tau:delays*tau
    x_hat_0 = [y_data(:,i); x_hat_0];
end
x_hat_0 = [x_data(2,delays*tau); x_hat_0];

t_cut = t(delays*tau:end); % Cut off first part to match model
x_data_cut = x_data(:,delays*tau:end);

% x_tilde_0 = U_hat'*x_hat_0; % Convert to tilde coordinate system
% x_tilde_hat = run_model(A_tilde,B_tilde,u_data,t_cut,x_tilde_0);
% x_hat = U_hat*x_tilde_hat; % Convert back to original coordinate system

x_hat = run_model(A,B,u_data,t_cut,x_hat_0);

plot(t_cut, x_data_cut)

hold on;
plot(t_cut, x_hat(1:2,:), '--');
hold off;

delays
tau
MSE_dmdc = mean(((x_data_cut(1,:) - x_hat(1,:)).^2)')'
%%
stop

%% Analytic system

A_c = [0, 1; -k/m, -b/m];
B_c = [0; 1/m];
C_c = [1 0];
D_c = 0;
sys_c = ss(A_c,B_c,C_c,D_c);
sys_d = c2d(sys_c, Ts);
[A_d,B_d,C_d,D_d] = ssdata(sys_d);

%% Plot and compare analytic to DMD model
x_hat_data = plot_model(A_d,B_d,u_data,t,x0);
plot(t, x_data-x_hat_data)
MSE_analytic = mean(((x_data-x_hat_data).^2)')'

%% Moving window
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
    Upsilon = U_dwork;
    
    % Inputs
    x       = x_data(k);
    u       = u_data(k);
    
    % Add input data to X2
    X2 = [X(:, 2:end), [x_dot; x]];
   
    % Based on DMD control example video by Steve Brunton
    XU = [X; Upsilon];
    G = X2*pinv(XU);
    A  = G(:,1:2);
    B  = G(:,end);
    
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
    U_dwork = [Upsilon(:, 2:end), u];
    
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
X = zeros(ny,1);
Upsilon = zeros(nu,1);
DX = [0; 0];

%% Local functions
function X_hat = plot_model(A,B,U_data,t,x0)
    N = max(size(t));   % Number of time steps 
    X_hat = zeros(length(x0),N); % Estimated X from model
    X_hat(:,1) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
    plot(t, X_hat);     % Estimated x
end

function X_hat = run_model(A,B,U_data,t,x0)
    N = max(size(t));   % Number of time steps 
    X_hat = zeros(length(x0),N); % Estimated X from model
    X_hat(:,1) = x0; % Initial conditions
    for index = 1:1:N-1
        X_hat(:,index+1) = A*X_hat(:,index) + B*U_data(index);
    end
end





