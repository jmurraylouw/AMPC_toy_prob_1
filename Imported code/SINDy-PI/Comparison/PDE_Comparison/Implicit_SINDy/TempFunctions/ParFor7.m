function f_x = ParFor7(in1)
%PARFOR7
%    F_X = PARFOR7(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    20-Sep-2019 09:35:36

u = in1(:,1);
ux = in1(:,4);
uxxx = in1(:,6);
f_x = ((u.*6.008914612862976e15+ux.*5.758587512457003e19-uxxx.*1.659241949896376e18-u.*ux.*4.555999807721701e19-u.*uxxx.*7.064407450052461e18+3.889391410438144e15).*-1.0)./(u.*3.455797482160128e15-1.115008480599081e17);
