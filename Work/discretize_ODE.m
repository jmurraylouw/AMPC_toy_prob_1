% Convert ODE to Difference Equation
% Can use different methods, like Tustin, Centred FD, Forward- or Backward Euler

disp("discretize_ODE")
disp("--------------")

syms x(t) th(t) u(t) % time dependant
syms s z 
syms X_z Th_z U_z % Z transform variables
syms p1 p2 p3 p4 % parameters (coef of terms)

% Define system as ODE
T = 0.01; % Sampling time
Dx = diff(x,t); % x_dot
D2x = diff(Dx, t); % x_dotdot

Dth = diff(th,t); % theta_dot
D2th = diff(Dth, t); % theta_dotdot

ODE_sys = D2x == ( p1*cos(th)*sin(th) + p2*Dth^2*sin(th) + p3*Dx + p4*u )/(1 - cos(th)^2); % Ordinary Diff Eq of system

% extract numerators and denomenators of both sides
[num_L, denom_L] = numden(lhs(ODE_sys)); 
[num_R, denom_R] = numden(rhs(ODE_sys));

% Flatten all fractions
ODE_sys = 0 == expand(num_L*denom_R - num_R*denom_L);

% Laplace transform
L_sys = laplace(ODE_sys);
pretty(L_sys)
stop
% Possible substitutions
Tustin = ((2/T)*(z - 1)/(z + 1));
Back_Euler = (1 - z^-1)/T;
Forward_Euler = (z - 1)/T;
Centered_FD = (z - z^-1)/(2*T); % Centered Finite Difference

% Convert to Z-transform
syms Dx(t) Dth(t)  % Resest symbol
old = [laplace(x), laplace(th),  laplace(u), diff(x,t), diff(th,t)];
new = [X_z,        Th_z,         U_z,        Dx,        Dth];
Z_sys = subs(L_sys, old, new); % Substitute with symbols to simplify
Z_sys = subs(Z_sys, s, Tustin); % Sub 's' for 'z' expression
Z_sys = simplifyFraction(Z_sys);

pretty(Z_sys)
stop

% extract numerators and denomenators of both sides
[num_L, denom_L] = numden(lhs(Z_sys)); 
[num_R, denom_R] = numden(rhs(Z_sys));

% Flatten all fractions
Z_sys = 0 == expand(num_L*denom_R - num_R*denom_L);

Z_sys = isolate(Z_sys, X_z*z);
Z_sys = expand(Z_sys);
Z_sys = collect(Z_sys, [X_z, F_z, z]);

pretty(Z_sys)

