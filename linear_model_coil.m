    clear; clc; close all;

%% Continuous-time system
A = [-433.3333, -5.0,       0.0;
	  47.61906, -0.0053571, 0.56057;
	  0.0,      25.0,       0.0];
B = [1.0;
     0.0;
     0.0];
C = [0.0, 0.0, 1.0];   % output: angle
D = 0.0;

SS = ss(A, B, C, D);

%% Discretize system
Ts = 1e-3;                 % sampling time
SimT = 10;                  % simulation duration
N = round(SimT/Ts);
t = (0:N-1)*Ts;
disableInput = 0;

SSd = c2d(SS, Ts);
Ad = SSd.A;
Bd = SSd.B;
Cd = SSd.C;
Dd = SSd.D;

%% Reference input
ref_type = 'sin';        % 'sin', 'square', 'saw'
A_ref = 0.1;                  % amplitude
f_ref = 1;               % Hz
DC = 0;


switch ref_type
    case 'sin'
        u_ref = A_ref * sin(2*pi*f_ref*t) + DC;
    case 'square'
        % Manual square wave [-A_ref, A_ref]
        u_ref = A_ref * sign(sin(2*pi*f_ref*t)) + DC;
    case 'saw'
        % Manual sawtooth from -A_ref to +A_ref
        u_ref = A_ref * (2*(mod(f_ref*t,1)) - 1) + DC;
    case 'constant'
        u_ref = 0*t + DC;  
end

%% Filter variables
tau = 0.01;
alpha = tau/(tau+Ts);
A_noise = 5e-3;
Limit = 100;

%% State-feedback controller design
% desired_poles = [0.8 0.5 0.6];   % discrete-time closed-loop poles
% K = place(Ad, Bd, desired_poles);

qCon = diag([5e0, 5e1, 5e4]);
rCon = 1e-3;

K = dlqr(Ad, Bd, qCon, rCon);

% Reference precompensation (to track input)
Nbar = 1 / (Cd * ((eye(size(Ad)) - (Ad - Bd*K)) \ Bd));



%% Observer design (Luenberger)
% observer_poles = [0.1 0.12 0.13]; % faster than closed-loop poles
% Ld = place(Ad', Cd', observer_poles)';

% Process noise (model uncertainty)
qObs = 0*eye(3);  

% Measurement noise (sensor noise)
rObs = A_noise^2;  % variance, not std dev

% Compute discrete-time Kalman gain
[Ld, P, ~] = dlqe(Ad, eye(3), Cd, qObs, rObs);

%% Initialize variables
x = zeros(3,N);          % true states: [theta; theta_dot]
xhat = zeros(3,N);       % observer initial state
y = zeros(1,N);          % measured output
yn = zeros(1,N);
yf = zeros(1,N);         % filtered output
u = zeros(1,N);          % applied control input




%% Simulation loop
for i = 1:N-1
    y(i) = Cd*x(:,i);
    yn(i) = y(i) + randn * A_noise;

    % Filter update
    yf(i+1) = alpha*yf(i) + (1-alpha)*yn(i);

    
    switch disableInput
        case 0
        u(i) = max(-Limit, min(Limit, -K*xhat(:,i) + Nbar*u_ref(i)));
        % u(i) =  -K*xhat(:,i) + Nbar*u_ref(i);
        case 1
        u(i) = 0;
    end

    % Plant update
    x(:,i+1) = Ad*x(:,i) + Bd*u(i);

    % Observer update
    xhat(:,i+1) = observer_step(xhat(:,i), u(i), yf(i), Ad, Bd, Cd, Ld);
    
end


%% Plot results
figure('Position',[645 0 895 880]);

subplot(4,1,1)
plot(t, u_ref, 'k:', 'LineWidth', 1.5); hold on;
plot(t, u, 'r', 'LineWidth', 1.5);
ylabel('Input');
xlabel('Time [s]');
legend('Reference', 'Applied u');
grid on;
title('Control Input');

subplot(4,1,2)
plot(t, yn(), 'r', 'LineWidth', 1.5); hold on;
plot(t, u_ref, 'm', 'LineWidth', 1.5); hold on;
plot(t, y(), 'g', 'LineWidth', 1.5); hold on;
plot(t, yf(), 'b', 'LineWidth', 1.5); hold on;
plot(t, xhat(3,:), 'k--', 'LineWidth', 1.5);
ylabel('Angle \theta [rad]');
xlabel('Time [s]');
legend( 'Measured', 'Reference', 'Actual',  'Measured (filtered)', 'Estimated');
grid on;
title('Pendulum Angle');

subplot(4,1,3)
plot(t, x(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, xhat(2,:), 'r:', 'LineWidth', 1.5);
ylabel('Angular Velocity \thetȧ [rad/s]');
xlabel('Time [s]');
legend('True', 'Estimated');
grid on;
title('Pendulum Angular Velocity');

subplot(4,1,4)
plot(t, x(1,:), 'g', 'LineWidth', 1.5); hold on;
plot(t, xhat(1,:),'k:', 'LineWidth', 1.5);
ylabel('Current I [A]');
xlabel('Time [s]')
legend('True', 'Estimated');
grid on;
title('Motor Current');


% Observer Function
 function xhat_next = observer_step(xhat, u, y, Ad, Bd, Cd, Ld)
     yhat = Cd * xhat;                        % estimated output (angle)
     xhat_next = Ad * xhat + Bd*u + Ld*(y - yhat); % Luenberger update
 end
