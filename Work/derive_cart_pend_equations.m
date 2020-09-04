syms m M l g d u
syms x(t) theta(t) xdot thetadot
y = sym('y', [4 1])

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

% Equation from Lagraungian, with regards to x
eq_x        = euler_lag(L, x,       u-d*dx,  t); % eq_x == 0
eq_theta    = euler_lag(L, theta,   0, t);

eq_x        = simplify(eq_x);
eq_theta    = simplify(eq_theta);

% Substitute symbols into derivatives
% Clear symbol connections
syms dx dtheta
syms ddx ddtheta

old = [diff(x,t,t), diff(theta,t,t)];
new = [ddx,    ddtheta];
eq_x = subs(eq_x, old, new);
eq_theta = subs(eq_theta, old, new);

old = [diff(x,t), diff(theta,t)];
new = [dx,    dtheta];
eq_x = subs(eq_x, old, new);
eq_theta = subs(eq_theta, old, new);

eq_x = isolate(eq_x==0, ddx);
eq_theta = isolate(eq_theta==0, ddtheta);

% Solve
eqns = [eq_x, eq_theta];
vars = [ddx, ddtheta];
solution = solve(eqns, vars);

solution.ddx
solution.ddtheta

ddx = simplifyFraction(solution.ddx);
ddtheta = simplifyFraction(solution.ddtheta);

% Substitute state variables with y
old = [x,    dx,   theta, dtheta];
new = [y(1), y(2), y(3),  y(4)];
ddx = subs(ddx, old, new);
ddtheta = subs(ddtheta, old, new);

ddx = simplifyFraction(solution.ddx);
ddtheta = simplifyFraction(solution.ddtheta);

ddtheta = simplifyFraction(subs(ddtheta, x, y(1)));
ddtheta = simplifyFraction(subs(ddtheta, theta, y(3)));

ddtheta
ddx

% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(thetadotdot);



