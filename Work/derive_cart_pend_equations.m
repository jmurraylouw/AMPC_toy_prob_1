clear all

syms m M l g d u
syms x(t) theta(t)
y = sym('y', [4 1]);
states = [x; theta];

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
PE_m = -m*g*z_m;

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
dstates  = [dx;  dtheta];
ddstates = [ddx; ddtheta];

eqns = [eq_x; eq_theta]; % Equations to solve with

% Substitute symbols into derivatives
old = [diff(states,t); diff(states,t,t)];
new = [dstates;        ddstates];
eqns = subs(eqns, old, new);

% Solve
solution = solve(eqns, ddstates);
ddstates = struct2cell(solution);
ddstates = [ddstates{:}]; % Convert to normal syms array from cell

% Simplify
ddstates = simplifyFraction(ddstates);

% Substitute state variables with y
old = [x;    dx;   theta; dtheta];
new = [y(1); y(2); y(3);  y(4)];
ddstates = subs(ddstates, old, new);

pretty(ddstates(1))
pretty(ddstates(2))

% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(thetadotdot);



