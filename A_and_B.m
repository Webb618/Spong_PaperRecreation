function [A,B] = A_and_B(q1_ref)

syms q1 q2 dq1 dq2 u
m1 = 1;   m2 = 1;
l1 = 1;   l2 = 2;
lc1 = 0.5; lc2 = 1;
I1 = 0.083; I2 = 0.33;
g  = 9.8;


d11 = m1*(lc1)^2 +m2*(l1^2+lc2^2+2*l1*lc2*cos(q2))+I1+I2;
d12 = m2*(lc2^2+l1*lc2*cos(q2)) +I2;
d21 = d12;
d22 = m2*lc2^2+I2;

    
h1 = -m2*l1*lc2*sin(q2)*dq2^2 - 2*m2*l1*lc2*sin(q2)*dq2*dq1;
h2 = m2*l1*lc2*sin(q2)*dq1^2;

    
phi1 = (m1*lc1+m2*l1)*g*cos(q1)+m2*lc2*g*cos(q1+q2);
phi2 = m2*lc2*g*cos(q1+q2);



 D = [d11 d12; d21 d22];
    H = [h1; h2];
    Phi = [phi1; phi2];

    ddq = D \ ( [0; u] - H - Phi );
    ddq1 = ddq(1);
    ddq2=ddq(2);


f1 = dq1;
f2 = dq2;
f3 = ddq1;
f4 = ddq2;
f= [f1;f2;f3;f4];
x= [q1;q2;dq1;dq2];

A_sym = jacobian(f, x);   % ∂f/∂x
B_sym = jacobian(f, u);   % ∂f/∂u

THETA4 =m1*lc1 +m1*l2;

x_equilibrium = [q1_ref;pi/2-q1_ref;0;0];
u_equilirbirum = 0;

A = double(subs(A_sym, {q1, q2,dq1,dq2, u}, {x_equilibrium(1),x_equilibrium(2),0,0,u_equilirbirum}));
B = double(subs(B_sym, {q1, q2,dq1,dq2, u}, {x_equilibrium(1),x_equilibrium(2),0,0,u_equilirbirum}));
end