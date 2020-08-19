    function [dx, y] = cartpend(t,x,u,M,m,L,d,g,varargin)
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

    % Output measurements
    y = zeros(2,1);
    y(1,1) = x(1);
    y(2,1) = x(3);
    
    % Derivatives
    dx = zeros(4,1);

    Sx = sin(x(3));
    Cx = cos(x(3));
    
    % Equations from derive_cartpend.m
    dx(1,1) = x(2);
    dx(2,1) = (L*m*sin(x(3))*x(4)^2 + u - d*x(2) + g*m*cos(x(3))*sin(x(3)))/(- m*cos(x(3))^2 + M + m);
    dx(3,1) = x(4);
    dx(4,1) = -(L*m*cos(x(3))*sin(x(3))*x(4)^2 + u*cos(x(3)) - d*x(2)*cos(x(3)) + g*m*sin(x(3)) + M*g*sin(x(3)))/(L*(- m*cos(x(3))^2 + M + m));
end




