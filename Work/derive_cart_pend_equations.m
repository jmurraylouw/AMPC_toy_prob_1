syms m M l g d u
syms x(t) theta(t)
y = sym('y', [4 1]);

% Rates
dx        = diff(x, t);
dtheta    = diff(theta, t);

% Drone body equations
KE_M = 0.5*M*(dx^2); % Kinetic energy
PE_M = 0; % Potential energy

% Payload equations
x_m = x + l*sin(theta);     % x of payload
z_m = - l*cos(theta);       % z of payload

KE_m = 0.5*m*( diff(x_m,t)^2 + diff(z_m,t)^2 );
PE_m = m*g*z_m;

% Lagraugian
L = (KE_M + KE_m) - (PE_M + PE_m);
L = simplify(L);

% Non-conservative Forces and Torques
Qx     = u-d*dx;
Qtheta = 0;

% Lagraungian equations
eq_x     = euler_lag(L, x,       Qx,  t); % eq_x == 0
eq_theta = euler_lag(L, theta,   Qtheta, t);

% Clear symbol connections
syms dx dtheta
syms ddx ddtheta

vars = [ddx; ddtheta]; % Variables to solve for
eqns = [eq_x; eq_theta]; % VEquations to solve with

% Substitute symbols into derivatives
old = [diff(x,t), diff(theta,t), diff(x,t,t), diff(theta,t,t)];
new = [dx,        dtheta,        ddx,         ddtheta];
eqns = subs(eqns, old, new);

% Solve
solution = solve(eqns, vars);

% Simplify
ddx = simplifyFraction(solution.ddx);
ddtheta = simplifyFraction(solution.ddtheta);

% Substitute state variables with y
old = [x,    dx,   theta, dtheta];
new = [y(1), y(2), y(3),  y(4)];
ddx = subs(ddx, old, new);
ddtheta = subs(ddtheta, old, new);

pretty(ddx)
pretty(ddtheta)

% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(thetadotdot);



