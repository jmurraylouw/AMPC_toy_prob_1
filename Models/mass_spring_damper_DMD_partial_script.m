%% Investigate DMD with partial state feedback for mass spring damper
% Key words:
% Time delay embedding
% Hankel matrix

%% Variables for simulation
x0 = [5; 0];
m = 1;
b = 0.1;
k = 5;

%% Read data
u_data  = out.u.Data';
x_data  = out.x.Data';
y_data  = [1 0]*out.x.Data';
% y_data  = out.y.Data';
t       = out.tout';

nx = size(x_data)*[1; 0];
ny = size(y_data)*[1; 0];
nu = size(u_data)*[1; 0];

%%
N       = max(size(x_data));
Ts      = t(2)-t(1);     % Sample time of data

%% Hankel
window = 5000; 
max_delay = 300
interval = 1

x = x_data(:,1:window+1)'; 
scale_data = zeros(nx*nx, max_delay);
S_data = zeros(nx, max_delay);
phase_data = zeros(nx, max_delay);

for delays = 2:interval:max_delay

H = zeros(delays, window+1);
for j = 1:1:delays
    H(j,:) = y_data(:,j:window+j);
end

[U,S,V] = svd(H, 'econ');
v = V(:,1:nx);

scale = v'*x; % value to scale v by to get x

scale_data(:,delays) = reshape(scale,[nx*nx,1]);
S_data(:,delays) = [S(1,1); S(2,2)];
phase_data(:,delays) = [lag1; lag2];

end

%%
delay_arr = 1:1:max_delay;
plot(delay_arr, S_data, '.', 'MarkerSize', 10); hold on
plot(delay_arr, scale_data, '.','MarkerSize', 10); hold off
legend('S1', 'S2', '1,1', '2,1', '1,2', '2,2')
%%


%% Plot V vs X
plot(t(1:end-delays+1), x_data(1,1:end-delays+1)); hold on;
plot(t(1:end-delays+1), x_data(2,1:end-delays+1)); 

plot(t(1:end-delays+1), a1.*V(:,1)', 'o'); 
plot(t(1:end-delays+1), a2.*V(:,2)', 'o'); hold off;
% !!! S(1,1).*V(:,2)' == x_data(2,1:end-delays+1)
% But only works well for delays = 5
legend('x', 'xdot', 'v1', 'v2')


%% Plot diag(s)

plot(diag(S), '.', 'MarkerSize', [20])
S(1,1)
S(2,2)
S(3,3)

%% Batch DMD - Full state feedback
X = x_data(:, 1:end-1);
X2 = x_data(:, 2:end);
U = u_data(:, 1:end-1);
% X2 = A*X + B*U
% X2 = [A, B]*[X; U]

AB = X2/[X; U];
A = AB(:, 1:nx);
B = AB(:, nx+1:end);

A_c = [0, 1; -k/m, -b/m];
B_c = [0; 1/m];
C_c = [1 0];
D_c = 0;
sys_c = ss(A_c,B_c,C_c,D_c);
sys_d = c2d(sys_c, Ts);
[A_d,B_d,C_d,D_d] = ssdata(sys_d);

x_hat_data = plot_model(A,B_d,u_data,t,x0);
plot(t, x_data-x_hat_data)

MSE_dmd = mean(((x_data-x_hat_data).^2)')'

function thetadeg = phase_shift(sig1,sig2)
    % Find phase shift angle in degrees between two signals
    % From: mathworks.com/matlabcentral/answers/439172-how-to-determine-the-phase-difference-phase-shift-between-two-signals
    sig1s = [mean(sig1); 2*std(sig1)];
    sig2s = [mean(sig2); 2*std(sig2)];
    sig_sum = sig1 + sig2;
    sig_sums = [mean(sig_sum); 2*std(sig_sum)];
    c_fcn = @(theta) sqrt(sig1s(2).^2 + sig2s(2).^2 + 2*sig1s(2).*sig2s(2).*cos(theta)) - sig_sums(2);
    theta = fzero(c_fcn, 1);
    thetadeg = theta*180/pi;
end