%% Derive system equations for 2D drone with suspended payload
% Two vertical forces at distance, r, from COM represent the rotor forces

% Define symbolic variables
syms M % Mass of drone body (at fulcrum)
syms m % Mass of swinging payload
syms I % Moment of inertia of drone body
syms l % Length of pendulum
syms r % Distance from each rotor force to COM of drone
syms g % Acceleration due to gravity (always negative)
syms cx % Damping coef. of drone through air in x direction (f = cx*xdot)
syms cz % Damping coef. of drone in z direction (f = cy*zdot)
syms ctheta % Damping coef. of drone in theta direction (f = ct*thetadot)
syms cbeta % Damping coef. of drone in theta direction (f = ct*thetadot)
syms t % Time

syms F1 % Rotor force on drone on left of COM
syms F2 % Rotor force on drone on right of COM

syms x(t) % x position of drone
syms z(t) % z position of drone
syms theta(t) % Pitch angle of drone (horisontal = 0 rad)
syms beta(t) % Suspended angle of payload cable (vertical down = 0 rad)
syms dx % dx/dt of drone 
syms dz % dz/dt of drone 
syms dtheta % dtheta/dt of drone
syms dbeta % dbeta/dt of payload cable

X = sym('X', [8 1]); % State vector [x; z; theta; beta; dx; dz; dtheta; dbeta]

% Rates
dx        = diff(x, t);
dz        = diff(z, t);
dtheta    = diff(theta, t);
dbeta     = diff(beta, t);

% Drone body equations
KE_M = 0.5*M*(dx^2 + dz^2) + 0.5*I*dtheta^2; % Kinetic energy of drone body (linear + rotational)
PE_M = -M*g*z; % Potential energy of drone body

% Payload equations
x_m = x + l*sin(beta); % x position of payload
z_m = z - l*cos(beta); % z position of payload

KE_m = 0.5*m*( diff(x_m,t)^2 + diff(z_m,t)^2 ); % Kinetic energy of payload
PE_m = -m*g*z_m; % Potential energy of payload

% Lagrangian
L = (KE_M + KE_m) - (PE_M + PE_m);
L = simplify(L);

% Non-conservative Forces 
Qx = -(F1 + F2)*sin(theta) - cx*dx; % ?? change cx air damping according to pitch
Qz = (F1 + F2)*cos(theta) - cz*dz;

% Non-conservative Torques
Qtheta = F2*r - F1*r - ctheta*dtheta; % Torques caused be rotor forces and air damping
Qbeta  = -cbeta*dbeta; % Torques caused air damping on rotation of cable

% Lagrangian equations (eq_x == 0)
eq_x     = euler_lag(L, x, Qx, t); 
eq_z     = euler_lag(L, z, Qz, t);
eq_theta = euler_lag(L, theta, Qtheta, t);
eq_beta  = euler_lag(L, beta, Qbeta, t);

% Isolate x"
ddx     = rhs(isolate(eq_x==0, diff(x,t,t))); % Still depends on z", beta"

% Substitute x"
eq_z        = subs(eq_z, diff(x,t,t), ddx);  % Contains z", beta"
eq_beta    = subs(eq_beta, diff(x,t,t), ddx);  % Contains z", beta"

% Isolate z"
zdotdot     = rhs(isolate(eq_z==0, diff(z,t,t))); % Still depends on beta"

% Substitute z"
eq_x        = subs(eq_x, diff(z,t,t), zdotdot);  % Contains x", beta"
eq_beta    = subs(eq_beta, diff(z,t,t), zdotdot);  % Contains only beta"

% Isolate beta"
betadotdot     = rhs(isolate(eq_beta==0, diff(beta,t,t))); % No more dotdot dependencies

% Substitute beta"
eq_x    = subs(eq_x, diff(beta,t,t), betadotdot);  % Contains only x"
eq_z    = subs(eq_z, diff(beta,t,t), betadotdot);  % Contains only z"

% Isolate x" and z"
ddx     = rhs(isolate(eq_x==0, diff(x,t,t))); % No more dotdot dependencies
zdotdot     = rhs(isolate(eq_z==0, diff(z,t,t))); % No more dotdot dependencies

% Simplify
ddx = simplifyFraction(ddx);
zdotdot = simplifyFraction(zdotdot);
betadotdot = simplifyFraction(betadotdot);

% Substitute vecolcity state variables with y
old = [diff(x,t),   diff(z,t),  diff(beta,t)];
new = [X(4),        X(5),       X(6)];
ddx     = subs(ddx, old, new);
zdotdot     = subs(zdotdot, old, new);
betadotdot = subs(betadotdot, old, new);

% Substitute other state variables seperately to avoid losing derivatives
old = [x,    z,    beta];
new = [X(1), X(2), X(3)];  
ddx     = subs(ddx, old, new);
zdotdot     = subs(zdotdot, old, new);
betadotdot = subs(betadotdot, old, new);

%% Display to copy
ddx
zdotdot
betadotdot

%% Display pretty equations
'xdotdot'
pretty(ddx)
'zdotdot'
pretty(zdotdot)
'betadotdot'
pretty(betadotdot)




% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(betadotdot);



