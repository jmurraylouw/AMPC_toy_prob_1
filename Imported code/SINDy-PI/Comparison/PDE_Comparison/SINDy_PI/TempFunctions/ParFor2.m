function f_x = ParFor2(in1)
%PARFOR2
%    F_X = PARFOR2(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Jul-2019 11:57:35

u = in1(:,1);
ux = in1(:,4);
uxx = in1(:,5);
uxxx = in1(:,6);
f_x = (u.*2.6068e4+ux.*7.4257e4+uxx.*2.6372e5+uxxx.*4.3251e5-u.*ux.*1.0e5+u.*uxx.*2.3836e4-u.*uxxx.*1.0264e6)./(u.*2.6078e4+2.6069e5);
