clear;
clear functions

%% HELPERS: EOM, ZeroDynamics, A_and_B, acrobot_dynamics, hybrid_control

%% EoM derivation
%Figure 1
Spong_EoMs = EOM();
fprintf('\n');
fprintf('Equation 1 = \n');
disp(Spong_EoMs(1))
fprintf('\n');
fprintf('Equation 2 = \n');
disp(Spong_EoMs(2))

%% Graph of zero dynamics
%Figure 2

y =ZeroDynamics();

%% PFL and LQR response

% parameters
p.m1 = 1;   p.m2 = 1;
p.l1 = 1;   p.l2 = 2;
p.lc1 = 0.5; p.lc2 = 1;
p.I1 = 0.083; p.I2 = 0.33;
p.g  = 9.8;



% swing-up target for q1
q1_ref = pi/2;

%% K matrix calculation for LQR 
[A,B] = A_and_B(q1_ref);

Q = [1 0 0 0; 
     0 1 0 0;
     0  0 1 0;
     0  0  0 1];
R = 1;

K = lqr(A,B,Q,R);
tspan = [0 10];

%% Figure 3 Recreation
%I can create the same figures as spong but using different gains!
x0 = [-pi/2; 0 ;0;0];
kp =20;
kd =6.5;
%                               time, states, controller, params

odefun = @(t,x) acrobot_dynamics(t, x, hybrid_control(x,p,q1_ref,kp,kd,K,t,1,false,false), p);
[t,X] = ode45(odefun, tspan, x0);
q1  = X(:,1);
q2  = X(:,2);
dq1 = X(:,3);
dq2 = X(:,4);
figure;
plot(t, q1,t, mod(q2,-2*pi) ); hold on;
yline(q1_ref, '--', 'Q_1 Reference');
yline(pi/2 -q1_ref, '--', 'Q_2 Reference');
title('PFL Response 1');

%% Figure 4 Recreation
tspan = [0 10];
x0 = [-pi/2; 0 ;0;0];
kp =15;
kd =6.5;
odefun = @(t,x) acrobot_dynamics(t, x, hybrid_control(x,p,q1_ref,kp,kd,K,t,1,false,false), p);
[t,X] = ode45(odefun, tspan, x0);
q1  = X(:,1);
q2  = X(:,2);
dq1 = X(:,3);
dq2 = X(:,4);
figure;
plot(t, q1,t, mod(q2,-2*pi) ); hold on;
yline(q1_ref, '--', 'Q_1 Reference');
yline(pi/2 -q1_ref, '--', 'Q_2 Reference');
title('PFL Response 2');

%% Figure 5 Recreation
tspan = [0 10];
x0 = [-pi/2; 0 ;0;0];
kp =20;
kd =6.5;
odefun = @(t,x) acrobot_dynamics(t, x, hybrid_control(x,p,q1_ref,kp,kd,K,t,1,false,true), p);
[t,X] = ode45(odefun, tspan, x0);
q1  = X(:,1);
q2  = X(:,2);
dq1 = X(:,3);
dq2 = X(:,4);
figure;
plot(t, q1,t, mod(q2,-2*pi) ); hold on;
yline(q1_ref, '--', 'Q_1 Reference');
yline(pi/2 -q1_ref, '--', 'Q_2 Reference');
ylabel('q1 (rad)');
title('PFL Response 1 with LQR');
clear functions
%% Figure 6 Recreation
x0 = [-pi/2-.4; -0.4 ;0;0];
kp =50;
kd =2.5;
odefun = @(t,x) acrobot_dynamics(t, x, hybrid_control(x,p,q1_ref,kp,kd,K,t,2, false, true), p);
tspan = [0 16];
[t,X] = ode45(odefun, tspan, x0);
q1  = X(:,1);
q2  = X(:,2);
dq1 = X(:,3);
dq2 = X(:,4);
figure;
plot(t, q1,t, mod(q2+pi,2*pi)-pi ); hold on;
yline(q1_ref, '--', 'Q_1 Reference');
yline(pi/2 -q1_ref, '--', 'Q_2 Reference');
title('PFL for the noncollocated case with LQR');

%% Figure 7 Recreation

E = zeros(length(t),1);

for k=1:length(t)
    
    d11 = p.m1*p.lc1^2 + p.m2*(p.l1^2 + p.lc2^2 + 2*p.l1*p.lc2*cos(q2(k))) + p.I1 + p.I2;
    d22 = p.m2*p.lc2^2 + p.I2;
    d12 = p.m2*(p.lc2^2 + p.l1*p.lc2*cos(q2(k))) + p.I2;

    % Kinetic energy
    T = 0.5*( d11*dq1(k)^2 + 2*d12*dq1(k)*dq2(k) + d22*dq2(k)^2 );

    % Potential energy
    V = (p.m1*p.lc1 + p.m2*p.l1)*p.g*sin(q1(k)) + p.m2*p.lc2*p.g*sin(q1(k) + q2(k));

    E(k) = T + V;
end

figure;
plot(t,E), grid on
ylabel('Energy')
title('Total energy during swing up')

