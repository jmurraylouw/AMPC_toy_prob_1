syms m M l g d u
syms x(t) theta(t) xdot thetadot
y = sym('y', [4 1])

xdot        = diff(x, t);
thetadot    = diff(theta, t);

% Drone body equations
KE_M = 0.5*M*(xdot^2); % Kinetic energy
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
eq_x        = euler_lag(L, x,       u-d*xdot,  t); % eq_x == 0
eq_theta    = euler_lag(L, theta,   0, t);

eq_x        = simplify(eq_x);
eq_theta    = simplify(eq_theta);

xdotdot     = rhs(isolate(eq_x==0, diff(x,t,t)));
eq_theta    = subs(eq_theta, diff(x,t,t), xdotdot);  % Substitute xdotdot in

thetadotdot = (simplify(rhs(isolate(eq_theta==0, diff(theta,t,t)))));
xdotdot     = (simplify(subs(xdotdot, diff(theta,t,t), thetadotdot)));

% Substitute state variables with y
old = [diff(x,t), diff(theta,t)];
new = [y(2),    y(4)];

xdotdot = subs(xdotdot, old, new);
thetadotdot = subs(thetadotdot, old, new);

xdotdot = simplifyFraction(subs(xdotdot, x, y(1)));
xdotdot = simplifyFraction(subs(xdotdot, theta, y(3)));

thetadotdot = simplifyFraction(subs(thetadotdot, x, y(1)));
thetadotdot = simplifyFraction(subs(thetadotdot, theta, y(3)));

thetadotdot
xdotdot

% dy(1,1) = y(2);
% dy(2,1) = double(xdotdot);
% dy(3,1) = y(4);
% dy(4,1) = double(thetadotdot);



