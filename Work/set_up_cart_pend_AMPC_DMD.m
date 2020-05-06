% %load('C:\Users\Murray\OneDrive - Stellenbosch University\Masters\AMPC_toy_prob_1\Data\mpc1_cart_pendulem.mat');
% sys_c = mpc1.Model.Plant; % Continuous system from mpc designed for cart pendulum
% sys_d = c2d(sys_c, mpc1.Ts, 'zoh') % Convert to discrete

%
%mpc1.Model.Plant = sys_d; % Set AMPC model to discrete
[A,B,C,D] = ssdata(mpc1.Model.Plant); % Get state space data
