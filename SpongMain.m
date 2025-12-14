clear;
clear functions

%% HELPERS: EOM, ZeroDynamics, A_and_B, acrobot_dynamics, hybrid_control

%% EoM derivation at the bottom

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
% For LQR we used plugged into the symbolic A and B matricies and we used
% the same Q and R as spong
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
%Some part of my dynamics or controller is slightly different from spong
%meaning I need to use different gain values
% I was not able to figure out which 
% For the non-collocated we take advantage of the zero dynamics and we can
% set our q1ref to a constant at our goal and the system will naturally swing up 
%Here Q2 swings around after q1 settles but does not make a full rotation
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

% Here we cause q2 to make a full rotation while q1 is in it target
% position
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
%Here we use LQR when we reach close to the top in order to stabalize at
%our equilibrium point. In order to best do this we want to start early as
%q2 is swinging up in order to best "catch it" with LQR. This stratagy
%works in both the collocated and non-collocated case
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
% Here we use LQR with the collocated case. Once again it best to
% switch to LQR as q2 is coming up toward our target equilibrium position,
% when its energy is high. With the collocated case we have torque on q2 in
% the direction of the velocity of q1 in order to pump energy into the
% system in each swing. We can not rely on the zero dynamics like in the
% non-collocated case. When both are close to and approaching equilibrium 
% we can activate LQR to balance the pendulum at our desired 
% equilibrium points. 

% Gains we found using the same system I used in the final project to find
% optimal gains for each angle. Here I found the gains that got q1 closest
% to its goal while dq1 was very small. 
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
%Here we show that we pump energy into the system after each swing. This
%shows that our collocated controller acts as we expect it to, in order to
%bring q1 up to equilibrium
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

%% EOM Derivation
%Figure 1
Spong_EoMs = EOM();
fprintf('\n');
fprintf('Equation 1 = \n');
disp(Spong_EoMs(1))
fprintf('\n');
fprintf('Equation 2 = \n');
disp(Spong_EoMs(2))



