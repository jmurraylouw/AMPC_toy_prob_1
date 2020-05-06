syms v(t) x(t) f(t) m b k T s z X_z F_z V_z n VAR

% Define system as ODE
Dx = diff(x,t);
D2x = diff(Dx, t);
ODE_sys = f == m*D2x + b*Dx + k*x; % Ordinary Diff Eq of system

% Laplace transform
L_sys = laplace(ODE_sys);

% Convert to Z-transform with Bilinear/Tustin substitution
Tustin = ((2/T)*(z - 1)/(z + 1));
syms Dx(t) % Resest symbol
old = [laplace(x),  laplace(f), diff(x,t)];
new = [X_z,         F_z,        Dx];
Z_sys = subs(L_sys, old, new);
Z_sys = subs(Z_sys, s, Tustin);
Z_sys = simplifyFraction(Z_sys);

[num_L, denom_L] = numden(lhs(Z_sys)); % extract left numerator and denomenator
[num_R, denom_R] = numden(rhs(Z_sys));

% Flatten all fractions
Z_sys = 0 == num_L*denom_R - num_R*denom_L;

Z_sys = isolate(Z_sys, X_z*z^2);
Z_sys = expand(Z_sys);
Z_sys = collect(Z_sys, [X_z, F_z, z]);

pretty(Z_sys)
Z_sys

