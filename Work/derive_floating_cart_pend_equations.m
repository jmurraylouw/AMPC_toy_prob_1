%% Derive system equations for 2D floating cart pendulum

% Define symbolic variables
syms M % Mass of drone body (at fulcrum)
syms m % Mass of swinging payload
syms L % Length of pendulum
syms g % Acceleration due to gravity (always negative)
syms cx % Damping coef. of drone in x direction (f = cx*xdot)
syms cz % Damping coef. of drone in z direction (f = cy*zdot)
syms ct % Damping coef. of drone in theta direction (f = ct*thetadot)
syms t % Time

syms fx % Input force on drone in x direction
syms fz % Input force on drone in y direction

syms x(t) % x position of drone
syms z(t) % z position of drone
syms theta(t) % Angle from horisontal of pendulum (down = 0)
syms xdot % dx/dt of drone 
syms zdot % dz/dt of drone 
syms thetadot % dtheta/dt of pendulum

y = sym('y', [6 1]); % State vector [x; z; theta; xdot; zdot; thetadot]

% Velocities
xdot        = diff(x, t);
zdot        = diff(z, t);
thetadot    = diff(theta, t);

% Drone body equations
KE_M = 0.5*M*(xdot^2 + zdot^2); % Kinetic energy of drone body
PE_M = -M*g*z; % Potential energy of drone body

% Payload equations
x_m = x + L*sin(theta); % x position of payload
z_m = z - L*cos(theta); % z position of payload

KE_m = 0.5*m*( diff(x_m,t)^2 + diff(z_m,t)^2 ); % Kinetic energy of payload
PE_m = -m*g*z_m; % Potential energy of payload

% Lagrangian
L = (KE_M + KE_m) - (PE_M + PE_m);
L = simplify(L);

% Lagrangian equation (eq_x == 0) with respect to x, z, theta
eq_x        = euler_lag(L, x,       fx - cx*xdot, t); % Q = Force and damping
eq_z        = euler_lag(L, z,       fz - cz*zdot, t); % Q = Force and damping
eq_theta    = euler_lag(L, theta,   -ct*thetadot, t); % Q = Damping

% % Simplify
% eq_x        = simplify(eq_x);
% eq_z        = simplify(eq_z);
% eq_theta    = simplify(eq_theta);

% Isolate x"
xdotdot     = rhs(isolate(eq_x==0, diff(x,t,t))); % Still depends on z", theta"

% Substitute x"
eq_z        = subs(eq_z, diff(x,t,t), xdotdot);  % Contains z", theta"
eq_theta    = subs(eq_theta, diff(x,t,t), xdotdot);  % Contains z", theta"

% Isolate z"
zdotdot     = rhs(isolate(eq_z==0, diff(z,t,t))); % Still depends on theta"

% Substitute z"
eq_x        = subs(eq_x, diff(z,t,t), zdotdot);  % Contains x", theta"
eq_theta    = subs(eq_theta, diff(z,t,t), zdotdot);  % Contains only theta"

% Isolate theta"
thetadotdot     = rhs(isolate(eq_theta==0, diff(theta,t,t))); % No more dotdot dependencies

% Substitute theta"
eq_x    = subs(eq_x, diff(theta,t,t), thetadotdot);  % Contains only x"
eq_z    = subs(eq_z, diff(theta,t,t), thetadotdot);  % Contains only z"

% Isolate x" and z"
xdotdot     = rhs(isolate(eq_x==0, diff(x,t,t))); % No more dotdot dependencies
zdotdot     = rhs(isolate(eq_z==0, diff(z,t,t))); % No more dotdot dependencies

% Simplify
xdotdot = simplifyFraction(xdotdot);
zdotdot = simplifyFraction(zdotdot);
thetadotdot = simplifyFraction(thetadotdot);

% Substitute vecolcity state variables with y
old = [diff(x,t),   diff(z,t),  diff(theta,t)];
new = [y(4),        y(5),       y(6)];
xdotdot     = subs(xdotdot, old, new);
zdotdot     = subs(zdotdot, old, new);
thetadotdot = subs(thetadotdot, old, new);

% Substitute other state variables seperately to avoid losing derivatives
old = [x,    z,    theta];
new = [y(1), y(2), y(3)];  
xdotdot     = subs(xdotdot, old, new);
zdotdot     = subs(zdotdot, old, new);
thetadotdot = subs(thetadotdot, old, new);

%% Display to copy
xdotdot
zdotdot
thetadotdot

%% Display pretty equations
'xdotdot'
pretty(xdotdot)
'zdotdot'
pretty(zdotdot)
'thetadotdot'
pretty(thetadotdot)




% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(thetadotdot);



