clear; clc; close all;

%% Plant

% linear system matrices
	A(1,1) = 0.0;
	A(1,2) = 0.6349206349206349;
	A(1,3) = 0.0;
	A(1,4) = 0.0;
	A(2,1) = 0.0;
	A(2,2) = -35.891284446945875;
	A(2,3) = -0.0250832251718952;
	A(2,4) = 0.6081841295788805;
	A(3,1) = 0;
	A(3,2) = -0.5747588132292484;
	A(3,3) = -0.010602760115950698;
	A(3,4) = 0.25708139156994725;
	A(4,1) = 0;
	A(4,2) = 0;
	A(4,3) = 59.1715976331361;
	A(4,4) = 0.0;
	B(1,1) = 0.0;
	B(2,1) = 6.989377037557411;
	B(3,1) = 0.11192706288504244;
	B(4,1) = 0.0;
	C(1,1) = 1;
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
Q_aug = diag([1e-6,1e-7,1e1,1e2]);

R = 3.16e-3;     

K_aug = dlqr(Ad, Bd, Q_aug, R);

Kx = K_aug(1:4);

Ob = obsv(Ad, Cd)
rand_Ob = rank(Ob)
det_ob = det(Ob)


%%  Kalman observer

Qn = 1e2;

Rn = 1e-6;

[kalmf, Ld, E] = kalman(sys_d, Qn, Rn, 0);
observer_poles = eig(Ad - Ld*Cd)
controller_poles = eig(Ad - Bd * Kx)