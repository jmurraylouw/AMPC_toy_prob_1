function Sindy_ODEs = Sindy_ODE_RHS(t,in2,u1)
%SINDY_ODE_RHS
%    SINDY_ODES = SINDY_ODE_RHS(T,IN2,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    21-Sep-2019 21:17:02

z1 = in2(:,1);
z3 = in2(:,3);
z4 = in2(:,4);
t2 = sin(z1);
t3 = t2.^2;
t4 = z1.*2.0;
t5 = sin(t4);
t6 = z3.^2;
Sindy_ODEs = [z3;z4;(t2.*4.8979e5-t5.*t6.*1.243e4)./(t3.*2.5307e4+2.4693e4);(t5.*-2.446e5+u1.*5.1553e4+t2.*t6.*4.9899e4)./(t3.*5.0209e4+4.9791e4)];
