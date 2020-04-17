f = out.f.Data';
x = out.x.Data';
x_dot = out.x_dot.Data';
%Theta_data = out.Theta.Data;
t = out.tout;
n = max(size(x));
Ts = 0.1;
Ts_c = 0.1;

X = [x; x_dot];
X2 = X(:,2:end);
X = X(:,1:end-1);
U = f(1:end-1);

XU = [X; U];
AB = X2*pinv(XU);

A = AB(:,1:2);
B = AB(:,end);
C = [0 1];
D = 0;
H_RLS = ss(A,B,C,D,Ts);
P_RLS = pole(H_RLS);
H_RLS = d2d(H_RLS, Ts_c, 'zoh'); % Resample to match controller's Ts
%[A,B,C,D] = ssdata(H_RLS)

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

%% Plot model vs measured
X_hat = zeros(2,n); % Estimated X from model
X_hat(:,1) = [0; 0]; % Initial conditions
for index = 1:1:n-1
    X_hat(:,index+1) = A*X_hat(:,index) + B*f(index);
end
x_hat = X_hat(1,:); % First row is x
plot(t, x, 'k'); hold on; % Measured x
pause
plot(t, x_hat, 'r');  % Estimated x
plot(t, x-x_hat);  % Plot the error of DMD results
hold off;
%% Continuous model deduced using modelling
Ac = [-0.5 -5; 1 0];
Bc = [1; 0];
Cc = [0 1];
Dc = 0;
Hc = ss(Ac, Bc, Cc, Dc);
Hd = c2d(Hc, Ts, 'zoh'); % Discrete model

P_d = pole(Hd);
%% Plot model predicted by RLS

Theta = Theta_data(n,:);
x_hat = zeros(1,n); % Inititial condition
for k=3:1:n-1
    x_hat(k) = Theta*[f(k-1); x(k-1); x(k-2)];
end
plot(t, x_hat, 'b'); hold on
plot(t, x-x_hat);  % Plot the error of RLS results
hold off;

legend('Measurements', 'DMD', 'DMD Error', 'RLS', 'RLS Error');
%%
[A,B,C,D] = ssdata(mpc1.Model.Plant)








