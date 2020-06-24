function dx = floating_pendulum_2D(t,x,u,params)
    % floating_pendulum_2D Models the continuous system equations of a 2D floating pendulem.
    %   x  = state vector [x; z; theta; x_dot; z_dot; theta_dot]
    %   dx = derivative of state vector
    %   u  = input vector [fx; fz]
    %   params = parameter vector
%     x0 = [0.5; 1; -0.2; -0.3; 0.2; 1]

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










