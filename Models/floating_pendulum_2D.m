function dx = floating_pendulum_2D(t,x,u,params)
    % floating_pendulum_2D Models the continuous system equations of a 2D floating pendulem.
    %   x  = state vector [x; z; theta; x_dot; z_dot; theta_dot]
    %   dx = derivative of state vector
    %   u  = input vector [fx; fz]
    %   params = parameter vector

    M  = params(1); % mass of base at pivot of pendulum
    m  = params(2); % mass of payload at end of pendulum
    L  = params(3); % length of pendulem rod
    g  = params(4); % acceleration due to gravity (always negative)
    cx = params(5); % Damping coef. of drone in x direction (f = cx*xdot)
    cz = params(6); % Damping coef. of drone in z direction (f = cy*zdot)
	ct = params(7); % Damping coef. of drone in theta direction (f = ct*thetadot)
    
    g = -abs(g); % Enforce negative gravity to avoid sign error
 
    Sin = sin(x(3)); % Compute once to save computations
    Cos = cos(x(3)); 
    
    dx = zeros(6,1);
    dx(1,1) = x(4); % x_dot
    dx(2,1) = x(5); % z_dot
    dx(3,1) = x(6); % theta_dot
    dx(4,1) = (L*fx + ct*x(6)*Cos - L*cx*x(4) + L^2*m*x(6)^2*Sin - L*g*m*Cos*Sin)/(L*(- m*Cos^2 + M + m)); % x_dotdot
    dx(5,1) = (fz + M*g - cz*x(5))/M; % z_dotdot
    dx(6,1) = -(M*ct*x(6) + ct*m*x(6) + L*fx*m*Cos - L*g*m^2*Sin + L^2*m^2*x(6)^2*Cos*Sin - L*M*g*m*Sin - L*cx*m*x(4)*Cos)/(L^2*m*(- m*Cos^2 + M + m)); % theta_dotdot
    
end










