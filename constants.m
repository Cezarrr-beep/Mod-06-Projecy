        clear; clc; close all;

%% Continuous-time system
   % linear system matrices
    A(1,1) = 1.0000055557074674;
	A(1,2) = 1.4077469779677672e-4;
	A(1,3) = -0.03724946297922129;
	A(2,1) = 0.33557140196434015;
	A(2,2) = 1.000023619915738;
	A(2,3) = -0.006249909895842499;
	A(3,1) = 9.371882804452956e-7;
	A(3,2) = 2.374718210590495e-5;
	A(3,3) = 0.9497320350096172;
	B(1,1) = -0.0012422683821324072;
	B(2,1) = -2.0843429230409523e-4;
	B(3,1) = -0.0016764350019361567;
	C(1,1) = 0.005305849739116081;
	C(1,2) = 0.03162315006534316;
	C(1,3) = -9.88197522084401e-5;
	D(1,1) = -3.2956355308312313e-6;

	SS = ss (A, B, C, D, 0.0);


%% Discretize system
Ts = 1e-3;

SSd = c2d(SS, Ts);
Ad = SSd.A;
Bd = SSd.B;
Cd = SSd.C;
Dd = SSd.D;

A_noise = 0.005;


%% State-feedback controller design
A_aug = [ Ad        zeros(3,1);
         -Cd        1 - 1e-10 ];

B_aug = [ Bd;
           0  ];

Q_aug = diag([1, 10, 1, 1]);

R = 1e2;

K_aug = dlqr(A_aug, B_aug, Q_aug, R);
Kx = K_aug(1:3);
Ki = K_aug(4); 
%% Observer design (Luenberger)
% Process noise (model uncertainty)
qObs = diag([1e-5, 1e-5, 1e-5]);
 
% Measurement noise (sensor noise)
rObs = A_noise^2;  % variance, not std dev

% Compute discrete-time Kalman gain
[Ld, P, ~] = dlqe(Ad, eye(3), Cd, qObs, rObs);





