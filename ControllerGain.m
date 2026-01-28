clear; clc; close all;

%% Plant

% linear system matrices
A(1,1) = 0.0;
	A(1,2) = 0.6349206349206349;
	A(1,3) = 0.0;
	A(1,4) = 0.0;
	A(2,1) = 0.0;
	A(2,2) = -36.02336976931349;
	A(2,3) = -0.03557416420269359;
	A(2,4) = 0.6672640046439018;
	A(3,1) = 0.0;
	A(3,2) = -0.5201529855513373;
	A(3,3) = -0.012403752903651196;
	A(3,4) = 0.23265698634395576;
	A(4,1) = 0.0;
	A(4,2) = 0.0;
	A(4,3) = 74.07407407407408;
	A(4,4) = 0.0;
	B(1,1) = 0.0;
	B(2,1) = 7.015098995995891;
	B(3,1) = 0.10129326351400329;
	B(4,1) = 0.0;
	C(1,1) = 1.0;
	C(1,2) = 0.0;
	C(1,3) = 0.0;
	C(1,4) = 1.0;
	D(1,1) = 0.0;
    


sys_c = ss(A, B, C, D);

%%  Discretization


Ts = 1e-3;

sys_d = c2d(sys_c, Ts, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

%% Augment with integrator
%  z[k+1] = z[k] + Ts*(r - y)


A_aug = [ Ad,           zeros(4,1);
         -Ts*Cd,        1 ];

B_aug = [ Bd;
           Ts ];

C_aug = [ Cd, 0 ];

Co = ctrb(Ad,Bd)
rank_Co = rank(Co)

n = size(Ad,1)

%%  Discrete LQR design

% State order: [x; x_dot; theta_dot; theta]
Q_aug = diag([1e-7,1e-7,1e-4,1e-3]);

R = 1e-6;     

K_aug = dlqr(Ad, Bd, Q_aug, R);

Kx = K_aug(1:4);

Ob = obsv(Ad, Cd)
rand_Ob = rank(Ob)
det_ob = det(Ob)


%%  Kalman observer

Qn = 1e2;

Rn = 1e-1;

[kalmfd, Ld, Ed] = kalman(sys_d, Qn, Rn, 0);
observer_poles = eig(Ad - Ld*Cd)
controller_poles = eig(Ad - Bd * Kx)
