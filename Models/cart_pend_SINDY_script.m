%% SINDY - Moving-Window of cart pendulum
% Estimate non-linear system model with moving window 
% of data from pendulum on a cart
% Full-state feedback

%% Generate Data
Beta = [10; 28; 8/3]; % Lorenz's parameters (chaotic)
n = 3;
x0=[-8; 8; 27];  % Initial condition
tspan=[0.01:0.01:50];
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) system(t,x,Beta),tspan,x0,options);

%% Compute Derivative
for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),Beta);
end

%% Build library and compute sparse regression
Theta = poolData(x,n,3);  % up to third order polynomials
lambda = 0.025;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n)
poolDataLIST({'x','y','z'},Xi,n,polyorder);

function dx = system(t,x)
    dx = [
            10*(x(2) - x(1));
            x(1)*(28 - x(3)) - x(2);
            x(1)*x(2) - (8/3)*x(3);
    ];
end


