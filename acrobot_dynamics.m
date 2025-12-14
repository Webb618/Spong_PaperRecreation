function dx = acrobot_dynamics(t, x, u, p)
    % Acrobot dynamics as defined by Spong

    q1  = x(1);  q2  = x(2);
    dq1 = x(3);  dq2 = x(4);
    
    m1 = p.m1;
    m2 =p.m2;
    lc1 = p.lc1;
    lc2 = p.lc2;
    l1 = p.l1;
    I1 =p.I1;
    I2 = p.I2;
    g=p.g;

    d11 = m1*(lc1)^2 +m2*(l1^2+lc2^2+2*l1*lc2*cos(q2))+I1+I2;
    d12 = m2*(lc2^2+l1*lc2*cos(q2)) +I2;
    d21 = d12;
    d22 = m2*lc2^2+I2;

    
    h1 = -m2*l1*lc2*sin(q2)*dq2^2 - 2*m2*l1*lc2*sin(q2)*dq2*dq1;
    h2 = m2*l1*lc2*sin(q2)*dq1^2;

    
    phi1 = (m1*lc1+m2*l1)*g*cos(q1)+m2*lc2*g*cos(q1+q2);
    phi2 = m2*lc2*g*cos(q1+q2);

    D1 =d21 -d22*d11/d12;
    H1 = h2 - h1*d22/d12;
    PHI1 = phi2 - phi1*d22/d12;
    
    D2 = d22-d21*d12/d11;
    H2 = h2 - h1 *d21/d11;
    PHI2 = phi2 - phi1*d21/d11;
    
    ddq2 = (u - PHI2 - H2)/D2;
    ddq1 = (u-PHI1 -H1)/D1;
   
    
    dx = [dq1; dq2; ddq1; ddq2];
end
