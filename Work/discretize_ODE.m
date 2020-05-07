% Convert ODE of Mass-Spring_Damper system to Difference Equation
% Can use different methods, like Tustin, Centred FD, Forward- or Backward Euler

syms v(t) x(t) f(t) m b k T s z X_z F_z V_z n VAR

% Define system as ODE
Dx = diff(x,t);
D2x = diff(Dx, t);
ODE_sys = f == m*D2x + b*Dx + k*x; % Ordinary Diff Eq of system

% Laplace transform
L_sys = laplace(ODE_sys);

% Convert to Z-transform
Tustin = ((2/T)*(z - 1)/(z + 1));
Back_Euler = (1 - z^-1)/T;
Forward_Euler = (z - 1)/T;
Centered_FD = (z - z^-1)/(2*T); % Centered Finite Difference

syms Dx(t) % Resest symbol
old = [laplace(x),  laplace(f), diff(x,t)];
new = [X_z,         F_z,        Dx];
Z_sys = subs(L_sys, old, new);
Z_sys = subs(Z_sys, s, Tustin); % Sub 's' for 'z' expression
Z_sys = simplifyFraction(Z_sys);

% extract numerators and denomenators of both sides
[num_L, denom_L] = numden(lhs(Z_sys)); 
[num_R, denom_R] = numden(rhs(Z_sys));

% Flatten all fractions
Z_sys = 0 == expand(num_L*denom_R - num_R*denom_L);

Z_sys = isolate(Z_sys, X_z*z^2);
Z_sys = expand(Z_sys);
Z_sys = collect(Z_sys, [X_z, F_z, z]);

pretty(Z_sys)

