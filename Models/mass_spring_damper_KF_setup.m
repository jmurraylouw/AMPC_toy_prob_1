clear all; close all;

nx = 2;
nu = 1;
ny = 1;

% System definition
m = 1;
b = 0.01;
k = 5;

A= [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;
Ts = 0.01;
sys_c = ss(A,B,C,D);  % Continuous system
sys_d = c2d(sys_c, Ts); % Discrete system
[F,G,H,D] = ssdata(sys_d);

sigma_a = 0.01; % Std dev of acceleration/force applied to model
Q = G*(sigma_a^2)*G';
R = 0.01;

% Initialise
x0 = [0; 0];
P0 = [0.1, 0.1; 0.1, 0.1];

disp('msd KF setup complete')





