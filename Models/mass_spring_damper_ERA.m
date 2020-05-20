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

sys_cont = ss(A_c,B_c,C_c,D_c);
[A,B,C,D] = ssdata(c2d(sys_cont, Ts))

%% Read data
u_data  = out.u.Data';

x_data  = out.x.Data';
% y_data  = C_c*x_data;
y_data  = out.y.Data';

t       = out.tout';
N       = length(t); % Number of data entries
Ts      = t(2)-t(1);
n = size(x_data)*[1; 0]; % Number of states
p = size(y_data)*[1; 0]; % Number of outputs
m = size(u_data)*[1; 0]; % Number of inputs

%% Steve Brunton OKID-ERA code
% r = rows in Hankel matrix
% s = columns in Hankel matrix (number of measurements)
% reduce = reduced rank
reduce = 2;
p = 200; 
[YY,M] = OKID(y_data, u_data, reduce, p);
mco = floor((length(YY)-1)/2); % mco = r = s
[A_til,B_til,C_til,D_til,Ur,Sr,Vr] = ERA(YY,mco,mco,m,p,reduce);
A_til
B_til
%% Get Transformation Matrix
% r = rows in Hankel matrix
% s = columns in Hankel matrix (number of measurements)
Vr_size = size(Vr); % Size of reduce V matrix from SVD
r = Vr_size(2); % reduced rank
s = Vr_size(1); % Number of data samples

W = zeros(r,s); % Controlability matrix with s matrix entries
W(:,1:p) = B;
prev_entry = B;
for i = m:m:(s-1)*m
    AnB =  A*prev_entry; % A^n*B
    W(:,(i+1):i+p) = AnB;
    prev_entry = AnB;
end

T = W*Vr*Sr^-0.5 % Transformation matrix to original state (x = T*x_til)


%% Run and plot model

x_til_hat = zeros(r,N);
y_hat = zeros(p,N);
for i = 2:1:N
    x_til_hat(:,i) = A_til*x_til_hat(:,i-1) + B_til*u_data(:,i-1);
    y_hat(:,i-1) = C_til*x_til_hat(:,i-1) + D_til*u_data(:,i-1) ;
end

x_hat = T*x_til_hat;
subplot(1,2,1), plot(t,x_data')
hold on;
subplot(1,2,1), plot(t,x_hat','--');
hold off;
subplot(1,2,2), plot(t,y_data')
% %% My ERA code starts here.
% %% Extended Controllabity matrix
% % r = rows in Hankel matrix
% % s = columns in Hankel matrix (number of measurements)
% r = 2;
% s = 1000;
% 
% W = zeros(n,s*m); % Controlability matrix with s matrix entries
% W(:,1:p) = B;
% prev_entry = B;
% for i = m:m:(s-1)*m
%     AnB =  A*prev_entry; % A^n*B
%     W(:,(i+1):i+p) = AnB;
%     prev_entry = AnB;
% end
% 
% %% Hankel matrix and SVD
% H0 = hankel(y_data(2:s*2)'); % Y(1) = y_data(2), Y = Markov parameter
% H0 = H0(1:r,1:s); % Make a r by s Hankel matrix
% 
% H1 = hankel(y_data(3:s*2)'); 
% H1 = H1(1:r,1:s); % Make a r by s Hankel matrix
% 
% [P,D,Q] = svd(H0, 'econ');
% 
% D_pow = D^-0.5; % Save computational time
% E_p = [eye(p); zeros(r-p, p)];
% E_m = [eye(m); zeros(m*s-m, m)];
% 
% %% Balanced model matrixes
% 
% A_til = D_pow*P'*H1*Q*D_pow % Balanced model system matrix A tilda
% B_til = D_pow*Q'*E_m % Balanced model input matrix B tilda
% C_til = E_p'*P*D_pow % Balanced model measurement matrix C tilda
% 
% T = W*Q*D_pow % Transformation matrix to original state (x = T*x_til)
% 
% %% run model
% 
% x_til_hat = zeros(r,s);
% y_hat = zeros(p,s);
% for i = 2:1:s
%     x_til_hat(:,i) = A_til*x_til_hat(:,i-1)+ B_til*u_data(:,i-1);
%     y_hat(:,i-1) = C_til*x_til_hat(:,i-1);
% end
% 
% % Plot model
% x_hat = 0.00125*T*x_til_hat;
% plot(x_hat','--');
% hold on;
% plot(x_data(:,1:s)')
% hold off;
