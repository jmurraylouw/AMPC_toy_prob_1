close all

m = 1;
M = 5;
L = 2;
g = -9.81;
d = 10;

tspan = 0:.1:20;
x0 = [0; 0; 0; 1];
% Theta = 0, downwards
% Theta_dot = 1, anti-clockwise
[t,x] = ode45(@(t,x)cartpend(x,m,M,L,g,d,0),tspan,x0);
plot(t,x);

%% Data from simulink
t = out.x.Time;
x = out.x.Data;
plot(t,x);

%%

t_pause = 0;
for k=1:length(t)
    drawcartpend_bw(x(k,:),m,M,L,t_pause);
end

hold off;
figure;
plot(t,x);
legend('x', 'xdot', 'theta', 'thetadot')


% function dy = pendcart(y,m,M,L,g,d,u)