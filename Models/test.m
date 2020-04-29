A = [1 2; 0 0]
B = [0; 0]
C = [1 0]
D=0
sys_c = ss(A,B,C,D)
sys_d = c2d(sys_c, 0.01, 'zoh')
