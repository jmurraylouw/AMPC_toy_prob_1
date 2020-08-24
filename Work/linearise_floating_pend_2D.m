%% linearise_floating_pend_2D
% determine A, B, C matrixes for linear model of floating pendulum in 2D

x0 = [0; 0; pi; 0];
u0 = 0;
f = @cartpend;

A = df_dx(f, x0, u0)
B = df_du(f, x0, u0)
% 
% h = 4*eps
% imag(cos(pi/3+h*1i))/h

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
    
    m = 1;
    M = 5;
    L = 2;
    g = -10;
    d = 1;
    
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