%  Linear Regression Estimate of full state based on N samples of partial
%  measurement.
% Derived from: Linear Systems with Sparse Inputs: Observability and Input Recovery
% http://www.vision.jhu.edu/assets/SefatiACC15.pdf
% x_N = A^N*pinv(O)*(Y_N - Gamma*U_N)
 
% Read simulation data
u_data = out.u.Data';
x_data = out.x.Data';
y_data = out.y.Data';
t = out.x.Time';

% Dimensions
[nx, n_time] = size(x_data);
[nu, n_time] = size(u_data);
[ny, n_time] = size(y_data);

% System parameters
m = 1;
b = 0.5;
k = 5;

% Define continuous system
A= [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = t(2)-t(1);
sys_c = ss(A,B,C,D);

% Convert to discrete System
sys_d = c2d(sys_c, Ts);
[A,B,C,D] = ssdata(sys_d)

O = obsv(A,C);
if rank(O)~=nx
    error("ERROR. System not observable")
end
pinv_O = pinv(O); % Calculate pseudoinverse once
Gamma = [];
for index = 1:1:nx
    new_col = [zeros(ny*(index-1), nu); ...
               D; ...
               O(1:(nx-index)*ny, :)*B; ];
    
    Gamma = [Gamma, new_col];
end

x_hat_data = zeros(nx,n_time); % Assign memory
for n = nx:1:n_time
    Y = y_data(:,(n-nx+1):n);
    Y = reshape(Y, nx*ny,1); % Form column vector of [[y0]; [y1]; [y2]; ...]
    U = u_data(:,(n-nx+1):n);
    U = reshape(U, nx*nu,1); % Form column vector
    
    x_hat = A^nx*pinv_O*(Y - Gamma*U);
    
    x_hat_data(:,n) = x_hat;
end

plot(t,x_data); hold on
pause
plot(t,x_hat_data,'x'); hold off
legend("x", "x_dot", "x_hat", "x_dot_hat")
