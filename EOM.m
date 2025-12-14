function Spong_EoMs = EOM()

clear all; format compact

fprintf('Acrobot Equations of Motion Derivation\n');


% Symbolic variables, angles: q1 and q2, masses m1 and m2, link lenghts l1
% and l2, COM positions l_c1 and l_c2, MOIs I1 and I2, gravity g, and torque
% tau.
syms q1 q2 m1 m2 l1 l2 l_c1 l_c2 I1 I2 g tau real

% 1a. Generalized Coordinates (GC's) and their derivatives
GC = [{q1}, {q2}];  % Using absolute angles HERE

% Time derivatives of q. dq1 and dq2 are angular velocities, d2q1 and d2q2
% are angular accelerations.
dq1 = fulldiff(q1, GC);
d2q1 = fulldiff(dq1, GC);
dq2 = fulldiff(q2, GC);
d2q2 = fulldiff(dq2, GC);

fprintf('Generalized Coordinates: q1 (shoulder), q2 (elbow)\n\n');



% Link 1 COM
x_c1 = l_c1 * sin(q1);
y_c1 = -l_c1 * cos(q1);  

% Link 2 COM
x_c2 = l1 * sin(q1) + l_c2 * sin(q1 + q2);
y_c2 = -l1 * cos(q1) - l_c2 * cos(q1 + q2);


%Velocities of COMs
dx_c1 = fulldiff(x_c1, GC);
dy_c1 = fulldiff(y_c1, GC);

dx_c2 = fulldiff(x_c2, GC);
dy_c2 = fulldiff(y_c2, GC);


% 2. Kinetic Energy


% Link 1: translational E + rotational E
T1 = (1/2) * m1 * (dx_c1^2 + dy_c1^2) + (1/2) * I1 * dq1^2;

% Link 2: translational E + rotational E
T2 = (1/2) * m2 * (dx_c2^2 + dy_c2^2) + (1/2) * I2 * (dq1 + dq2)^2;

% Total KE=T1+T2
T = simplify(T1 + T2);

fprintf('Kinetic Energy T = \n');
pretty(T)
fprintf('\n');

% 3. Potential Energy

%PE=Mgh
V1 = m1 * g * y_c1;
V2 = m2 * g * y_c2;

V = simplify(V1 + V2);

fprintf('Potential Energy V = \n');
pretty(V)
fprintf('\n');

% 4. Lagrangian L=T-V
L = T - V;
L = simplify(L);

%EOMs

% Equation 1: d/dt(dL/dq1') - dL/dq1 = 0 (no torque applied at at joint 1)
eq1 = fulldiff(diff(L, dq1), GC) - diff(L, q1);
eq1 = simplify(eq1);

fprintf('Equation 1 (Joint 1 - unactuated):\n');
disp(eq1);
fprintf('\n');

% Equation 2: d/dt(dL/dq2') - dL/dq2 = tau (torque at joint 2)
eq2 = fulldiff(diff(L, dq2), GC) - diff(L, q2);
eq2 = simplify(eq2);

fprintf('Equation 2 (Joint 2 - actuated):\n');
disp(eq2);
fprintf('\n');
Spong_EoMs = [eq1;eq2];

% 6. Xi: Non-conservative terms
Xi1 = 0;      % No applied torque at joint 1
Xi2 = tau;    % Torque at joint 2

fprintf('Non-conservative forces:\n');
fprintf('  Xi1 = 0 (no actuator at shoulder)\n');
fprintf('  Xi2 = tau (elbow torque)\n\n');

% 7. Get Mass Matrix D and Tau

% Mass matrix D
D = sym(zeros(2,2));
D(1,1) = diff(eq1 - Xi1, d2q1);
D(1,2) = diff(eq1 - Xi1, d2q2);
D(2,1) = diff(eq2 - Xi2, d2q1);
D(2,2) = diff(eq2 - Xi2, d2q2);

D = simplify(D);

fprintf('Mass Matrix D =\n');
disp(D);
fprintf('\n');

% Check symmetry of Mass Matrix D
fprintf('Checking symmetry: M(1,2) - M(2,1) = \n');
disp(simplify(D(1,2) - D(2,1)));
fprintf('(Should be zero for symmetric M)\n\n');

% Finding Tau 
Tau = sym(zeros(2,1));
Tau(1,1) = -(eq1 - Xi1) + D(1,1)*d2q1 + D(1,2)*d2q2;
Tau(2,1) = -(eq2 - Xi2) + D(2,1)*d2q1 + D(2,2)*d2q2;

Tau = simplify(Tau);

fprintf('Tau = [RHS of M*d2q = Tau] =\n');
disp(Tau);
fprintf('\n');

%Spong notation
fprintf('Converting to Spong (1995) notation...\n\n');

% Spong uses: d11, d12, d22, h1, h2, phi1, phi2
% Equation form: [d11 d12; d12 d22]*[d2q1; d2q2] + [h1; h2] + [phi1; phi2] = [0; tau]

% Extract mass matrix elements
d11 = D(1,1);
d12 = D(1,2);
d21 = D(2,1);
d22 = D(2,2);

fprintf('Mass Matrix Elements:\n');
fprintf('d11 = '); disp(d11);
fprintf('d12 = d21 = '); disp(d12);
fprintf('d22 = '); disp(d22);
fprintf('\n');

% Finding h terms (coriolis/centrifugal) and phi 

% Gravity terms (terms with g)
phi1 = -simplify(subs(Tau(1), [dq1, dq2], [0, 0]));
phi2 = -simplify(subs(Tau(2), [dq1, dq2], [0, 0]));

fprintf('Gravity Terms phi1 and phi2:\n');
fprintf('phi1 = '); disp(phi1);
fprintf('phi2 = '); disp(phi2);
fprintf('\n');

% Coriolis/centrifugal terms (remaining terms)
h1 = -Tau(1) - phi1;
h2 = -Tau(2) - phi2;

h1 = simplify(h1);
h2 = simplify(h2);

fprintf('Coriolis/Centrifugal Terms h1 and h2:\n');
fprintf('h1 = '); disp(h1);
fprintf('h2 = '); disp(h2);
fprintf('\n');


%Solving for accelerations

eqns = [eq1 - Xi1, eq2 - Xi2];
S = solve(eqns, [d2q1, d2q2]);

fprintf('Accelerations:\n');
fprintf('d2q1 = '); disp(S.d2q1);
fprintf('\nd2q2 = '); disp(S.d2q2);

end