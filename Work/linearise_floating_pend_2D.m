function [A_lin, B_lin] = linearise_floating_pend_2D(f,x0,u0)
%% linearise_floating_pend_2D
% determine A, B, C matrixes for linear model of floating pendulum in 2D

% x0 = [0; 0; 0; 0; 0; 0];
% u0 = [0; 6*9.81];
% f = @floating_pendulum_2D; % Non-linear model to linearise
f(x0, u0) % Check if f(x0, u0) = =0

A_lin = df_dx(f, x0, u0); % Linearised matrix
B_lin = df_du(f, x0, u0);
% 
% h = 4*eps
% imag(cos(pi/3+h*1i))/h

end

function A = df_dx(f, x0, u0)
    % Partial derivative of f with respect to x
    % JACCSD Jacobian through complex step differentiation
    % By Yi Cao at Cranfield University, 02/01/2008

    f_xu = f(x0,u0);
    n = numel(x0);
    m = numel(f_xu);
    A = zeros(m,n);
    h = n*eps*10; % Very small value

    for k = 1:n
        x1 = x0;
        x1(k) = x1(k)+ h*1i;
        A(:,k) = imag(f(x1,u0))/h;
    end

end

function B = df_du(f, x0, u0)
    % Partial derivative of f with respect to x
    % JACCSD Jacobian through complex step differentiation
    % By Yi Cao at Cranfield University, 02/01/2008

    f_xu = f(x0,u0);
    l = numel(u0);
    m = numel(f_xu);
    B = zeros(m,l);
    h = l*eps; % Very small value

    for k = 1:l
        u1 = u0;
        u1(k) = u1(k) + h*1i; % Add small imaginary part
        B(:,k) = imag(f(x0,u1))/h;
    end

end

function dx = floating_pendulum_2D(x,u)
    % floating_pendulum_2D Models the continuous system equations of a 2D floating pendulem.
    %   x  = state vector [x; z; theta; x_dot; z_dot; theta_dot]
    %   dx = derivative of state vector
    %   u  = input vector [fx; fz]
    %   params = parameter vector
%     x0 = [0.5; 6; -0.2; -0.8; -1; 1]

    M  = 4; % mass of base at pivot of pendulum
    m  = 2; % mass of payload at end of pendulum
    L  = 0.8; % length of pendulem rod
    g  = -9.81; % acceleration due to gravity (always negative)
    cx = 0.01; % Damping coef. of drone in x direction (f = cx*xdot)
    cz = 0.1; % Damping coef. of drone in z direction (f = cy*zdot)
	ct = 0.01; % Damping coef. of drone in theta direction (f = ct*thetadot)
    
    g = -abs(g); % Enforce negative gravity to avoid sign error
 
    Sin = sin(x(3)); % Compute once to save computations
    Cos = cos(x(3)); 
    
    fx = u(1);
    fz = u(2);
    
    dx = zeros(6,1);
    dx(1,1) = x(4); % x_dot
    dx(2,1) = x(5); % z_dot
    dx(3,1) = x(6); % theta_dot
    dx(4,1) = (M*ct*x(6)*Cos + ct*m*x(6)*Cos + L*M*fx*Cos^2 + L*M*fx*Sin^2 + L*fx*m*Cos^2 - L*M*cx*x(4)*Cos^2 - L*M*cx*x(4)*Sin^2 - L*cx*m*x(4)*Cos^2 + L*fz*m*Cos*Sin + L^2*M*m*x(6)^2*Sin^3 - L*cz*m*x(5)*Cos*Sin + L^2*M*m*x(6)^2*Cos^2*Sin)/((M + m)*(L*M*Cos^2 + L*M*Sin^2)); % x_dotdot
    dx(5,1) = (L*M^2*g*Sin^2 + M*ct*x(6)*Sin + ct*m*x(6)*Sin + L*M*fz*Cos^2 + L*M*fz*Sin^2 + L*fz*m*Sin^2 + L*M^2*g*Cos^2 + L*M*g*m*Cos^2 - L*M*cz*x(5)*Cos^2 + L*M*g*m*Sin^2 - L*M*cz*x(5)*Sin^2 - L*cz*m*x(5)*Sin^2 + L*fx*m*Cos*Sin - L^2*M*m*x(6)^2*Cos^3 - L*cx*m*x(4)*Cos*Sin - L^2*M*m*x(6)^2*Cos*Sin^2)/((M + m)*(L*M*Cos^2 + L*M*Sin^2)); % z_dotdot
    dx(6,1) = -(M*ct*x(6) + ct*m*x(6) + L*fx*m*Cos + L*fz*m*Sin - L*cx*m*x(4)*Cos - L*cz*m*x(5)*Sin)/(m*(M*L^2*Cos^2 + M*L^2*Sin^2)); % theta_dotdot
    
end

% For testing
function dx = cartpend(x,u)
    
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
    
    m = 2;
    M = 4;
    L = 1;
    g = -9.81;
    d = 5;
        
    % Derivatives
    dx = zeros(4,1);

    Sx = sin(x(3));
    Cx = cos(x(3));
    D = m*L*L*(M+m*(1-Cx^2));

    % Equations from derive_cartpend.m
    dx(1,1) = x(2);
    dx(2,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(2))) + m*L*L*(1/D)*u;
    dx(3,1) = x(4);
    dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(2))) - m*L*Cx*(1/D)*u; % +.01*randn;

end