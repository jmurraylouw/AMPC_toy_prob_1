function f_x = ParFor4(in1)
%PARFOR4
%    F_X = PARFOR4(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    31-Mar-2020 16:50:24

r = in1(:,1);
s = in1(:,13);
u = in1(:,19);
ut = in1(:,20);
ux = in1(:,21);
uxx = in1(:,22);
uy = in1(:,23);
uyy = in1(:,24);
z = in1(:,7);
t2 = u.^2;
f_x = r.*-4.845016406907234+s.*6.96612441301113+ut.*5.457984690783633e-1-ux.*1.530748504418966e-2+uxx.*2.135644249305187-uy.*3.570878181824355e-2+uyy.*5.02823932467436-z.*1.323504040227272e2+r.*u.*6.191151455437648-s.*u.*9.105578580027213-t2.*u.*1.077295023277402e2+u.*ut.*5.65024831920482e-2+u.*ux.*1.830355084143775e-2-u.*uxx.*4.132096419489244+u.*uy.*4.442896435307375e-2-u.*uyy.*6.867807677976089+u.*z.*1.83605552916415e2+4.028421339788474e1;
