%% SINDY - Moving-Window of cart pendulum
% Estimate non-linear system model with moving window 
% of data from pendulum on a cart
% Full-state feedback

% Adaptation of code by Steve Brunton

clear all;
%% Generate Data
x0 = [0; 0; 0.1; -4];  % Initial condition
n = length(x0);  % Number of states
tspan = [0.01:0.01:20];
u = square(tspan);

% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x] = ode45(@(t,x) cartpend(t,x,square(t)), tspan, x0);

%% Compute Derivative
dx = zeros(size(x));
for i=1:length(x)
    dx(i,:) = cartpend(0,x(i,:),square(t(i)));
end

%% Build library and compute sparse regression
polyorder = 3;

Sx = sin(x(:,3));
Cx = cos(x(:,3));
D = 1./(1-Cx.^2);

x_extended = [x, u', Sx, Cx, D]; % Add other terms to use in Theta
% x_extended = x;
Theta = poolData(x_extended,polyorder);  % up to third order polynomials

lambda = 0.015;      % lambda is our sparsification knob.
% Xi = sparsifyDynamics(Theta,dx,lambda,n)
Xi(:,1) = lasso(Theta,dx(:,1),'Lambda',lambda);
Xi(:,2) = lasso(Theta,dx(:,2),'Lambda',lambda);
Xi(:,3) = lasso(Theta,dx(:,3),'Lambda',lambda);
Xi(:,4) = lasso(Theta,dx(:,4),'Lambda',lambda);

figure
subplot(2,2,1), bar(Xi(:,1))
subplot(2,2,2), bar(Xi(:,2))
subplot(2,2,3), bar(Xi(:,3))
subplot(2,2,4), bar(Xi(:,4))

% Xi = lasso(Theta, dx(:,1), 'Lambda', lambda)
% subplot(1,2,2), bar(Xi)
% size(Xi)
x_names = {'x.','xdot.','th.','thdot.','Sx.','Cx.','D.'};
% visualize_Xi(x_names, Xi, polyorder) % Visualise terms of Xi

%% Run model

tspan = [0.01:0.01:40];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
u = square(tspan);
[t,x] = ode45(@(t,x) cartpend(t,x,square(t)), tspan, x0);
[t,x_hat] = ode45(@(t,x) sindy_system(t,x,square(t),Xi,polyorder), tspan, x0);

figure(2)
subplot(2,1,1), plot(t,x); hold on
subplot(2,1,1), plot(t,x_hat,'--')
subplot(2,1,2), plot(t,x-x_hat)

lambda
ave_error = mean(mean(abs(x-x_hat)))
%%

function dx = system(t,x)
    dx = [
            10*x(2) - 10*x(1);
            28*x(1) - x(3) - 2/(1 - sin(x(2)).^2);
            3 + x(1) - 3*x(3)*x(2);
    ];
end

function dx = cartpend(t,x,u)
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

    m = 1;
    M=1;
    L=1;
    g=-10;
    d=0.2;

    dx = zeros(4,1);

    Sx = sin(x(3));
    Cx = cos(x(3));
    D = m*L*L*(M+m*(1-Cx^2));
    
    dx(1,1) = x(2);
    dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
    dx(3,1) = x(4);
    dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;
end

function dx = sindy_system(t,x,u, Xi, polyorder)
    % Run model generated by SINDY.
    % x is a column vector
    % dx is a column vecotr
    x = x'; % Turn into row vector
    Sx = sin(x(:,3));
    Cx = cos(x(:,3));
    D = 1./(1-Cx.^2);

    x_extended = [x, u, Sx, Cx, D]; % Add other terms to use in Theta
    % x_extended = x;
    Theta = poolData(x_extended,polyorder);  % up to third order polynomials

    dx = Theta*Xi;
    dx = dx'; % Turn into column vector
end