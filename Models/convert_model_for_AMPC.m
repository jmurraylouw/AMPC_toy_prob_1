f = out.f.Data';
x = out.x.Data';
x_dot = out.x_dot.Data';
Theta_data = out.Theta.Data;
t = out.tout;
n = max(size(x));


% Mass spring damper parameters
m = 1;
b = 0.5;
k = 5;
Ts = 0.1;

% % Check hand calculated model
% a = b/(2*m);
% omega = sqrt(k/m - a^2);
% K = 1/(m*omega);
% 
% 2*exp(-a*Ts*cos(omega*Ts));
% exp(-2*a*Ts);

% Define model as State Space
A = [-b/m -k/m; 1 0];
B = [1/m; 0];
C = [0 1];
D = [0];
Hssc = ss(A,B,C,D); % Define continuous state space system
Hssd = c2d(Hssc, Ts, 'zoh'); % Convert to discrete system
[A,B,C,D] = ssdata(Hssd)
% Extract nx, ny, nu
size_B = size(B);
nx = size_B(1);
nu = size_B(2);
size_C = size(C);
ny = size_C(1);

%Create nominal point at all 0, becasue linear model.
X = zeros(nx,1);
Y = zeros(ny,1);
U = zeros(nu,1);
DX = [0; 0];

%% Plot model predicted by RLS
plot(t, x, 'k'); hold on
increments = 1/0.01
for index = 1:increments:n
    Theta = Theta_data(index,:);
    x_hat = zeros(1,n); % Inititial condition
    for k=3:1:n-1
        x_hat(k) = Theta*[f(k-1); x_hat(k-1); x_hat(k-2)];
    end
    
    plot(t, x_hat);
    pause
      % Plot the error of RLS results
end
plot(t, x-x_hat);
%% Convert RLS model to A,B,C,D

A = [Theta(2) Theta(3); 1 0];
B = [Theta(1); 0];
C = [1 0];
D = [0];

X_hat = zeros(2,n); % Estimated X from model
X_hat(:,1) = [0; 0]; % Initial conditions
for index = 1:1:n-1
    X_hat(:,index+1) = A*X_hat(:,index) + B*f(index);
end
x_hat = X_hat(1,:); % First row is x
plot(t, x_hat, 'r');  % Estimated x
plot(t, x-x_hat);  % Plot the error of DMD results
hold off;

legend('Measurements', 'Theta', 'Theta error', 'A,B', 'A,B error');



