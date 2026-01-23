clear; clc; close all;

%% Plant

A = [ -0.018122374991498984,   0.14122799209516126,  -38.20972164751147;
      335.5704697986577,       0.0,                  0.0;
     -0.003125793391744057,   0.02435936374936598,  -51.56395461829034 ];

B = [ -0.040296682154882524;
       0.0;
      -0.054380304286706566 ];

C = [0 1 0];
D = 0;

sys_c = ss(A, B, C, D);

%%  Discretization


Ts = 1e-3;   % 1000 Hz

sys_d = c2d(sys_c, Ts, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

%% Augment with integrator
%  z[k+1] = z[k] + Ts*(r - y)


A_aug = [ Ad,           zeros(3,1);
         -Ts*Cd,        1 ];

B_aug = [ Bd;
           Ts ];

C_aug = [ Cd, 0 ];

%%  Discrete LQR design

% State order: [x; x_dot; theta; integral]
Q_aug = diag([1, 10, 100, 500]);

R = 1;     

K_aug = dlqr(A_aug, B_aug, Q_aug, R);

Kx = K_aug(1:3);
Ki = K_aug(4);

disp('State feedback gain Kx = ');
disp(Kx)

disp('Integral gain Ki = ');
disp(Ki)

%%  Kalman observer

% Process noise
Qn = diag([1e-4, 1e-3, 1e-3]);

% Measurement noise
Rn = 1e-3;

[Ld, P, E] = dlqe(Ad, eye(3), Cd, Qn, Rn);

disp('Observer gain Ld = ');
disp(Ld)

%%  Closed-loop check

Acl = [Ad - Bd*Kx,    Bd*Ki;
       -Ts*Cd,       1 ];

eig_cl = eig(Acl);

disp('Closed-loop eigenvalues:');

Nbar = 1 / (Cd * ((eye(size(Ad)) - Ad + Bd*Kx) \ Bd));
