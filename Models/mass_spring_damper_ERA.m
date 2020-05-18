% Apply ERA to mass spring damper.
% Based on An Eigensystem Realization Algorithm
% for Modal Parameter Identification and Model Reduction
% by Juang and Pappa, 1985
%% Variables for simulation
x0 = [0; 0];
m = 1;
b = 0.1;
k = 5;

A_c = [0, 1; -k/m, -b/m];
B_c = [0; 1/m];
C_c = [1 0];
D_c = 0;

%% Read data
u_data  = out.u.Data';

x_data  = out.x.Data';
y_data  = C_c*x_data;
% y_data  = out.y.Data';
t       = out.tout';
Ts      = t(2)-t(1);
n = size(x_data)*[1; 0]; % Number of states
p = size(y_data)*[1; 0]; % Number of outputs
m = size(u_data)*[1; 0]; % Number of inputs

sys_cont = ss(A_c,B_c,C_c,D_c);
[A,B,C,D] = ssdata(c2d(sys_cont, Ts))

%% Extended Controllabity matrix
% r = rows in Hankel matrix
% s = columns in Hankel matrix (number of measurements)
r = 2;
s = 1000;

W = zeros(n,s*m); % Controlability matrix with s matrix entries
W(:,1:p) = B;
prev_entry = B;
for i = m:m:(s-1)*m
    AnB =  A*prev_entry; % A^n*B
    W(:,(i+1):i+p) = AnB;
    prev_entry = AnB;
end

%% Hankel matrix and SVD
H0 = hankel(y_data(1:s*2)'); 
H0 = H0(1:r,1:s); % Make a r by s Hankel matrix

H1 = hankel(y_data(2:s*2)'); 
H1 = H1(1:r,1:s); % Make a r by s Hankel matrix

[P,D,Q] = svd(H0, 'econ');

D_pow = D^-0.5; % Save computational time
E_p = [eye(p); zeros(n-p, p)];
E_m = [eye(m); zeros(m*s-m, m)];

%% Balanced model matrixes

A_til = D_pow*P'*H1*Q*D_pow % Balanced model system matrix A tilda
B_til = D_pow*Q'*E_m % Balanced model input matrix B tilda
C_til = E_p'*P*D_pow % Balanced model measurement matrix C tilda

T = W*Q*D_pow % Transformation matrix to original state (x = T*x_til)

%% run model

x_til_hat = zeros(n,s);
y_hat = zeros(p,s);
for i = 2:1:s
    x_til_hat(:,i) = A_til*x_til_hat(:,i-1)+ B_til*u_data(:,i-1);
    y_hat(:,i-1) = C_til*x_til_hat(:,i-1);
end
x_hat = T*x_til_hat;
plot(x_hat','--');
hold on;
plot(x_data(:,1:s)')
hold off;


