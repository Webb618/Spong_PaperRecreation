function u = hybrid_control(x, p, q1_ref, kp, kd,K,t, SIGMA, resetLatch, LQR)
    q1  = x(1);  q2  = x(2);
    dq1 = x(3);  dq2 = x(4);
    persistent latched
    if isempty(latched) || (nargin >= 20 && resetLatch)
        latched = false;
    end
    q1_wrapped = mod((q1 )+pi, 2*pi) -pi;
    q2_wrapped = mod((q2 )+pi, 2*pi) -pi;
    x_ref = [q1_ref;pi/2 - q1_ref; 0; 0];
    dx1_wrapped = q1_wrapped-x_ref(1);
    dx2_wrapped = q2_wrapped-x_ref(2);
    dx_wrapped = [dx1_wrapped;dx2_wrapped;x(3);x(4)];
    
    LQR_region = (abs(dx1_wrapped )< 0.6) && ...
                    (abs(dx2_wrapped) < 0.8)  && abs(dq2)< 2.0 && abs(dq1) < 1.5; 
    if ~latched
        if LQR_region && LQR
            latched = true;
        end
    end
    if latched
        
        u = -K * dx_wrapped;

        return
    end

     
        % model terms
    m1 = p.m1;
    m2 =p.m2;
    lc1 = p.lc1;
    lc2 = p.lc2;
    l1 = p.l1;
    I1 =p.I1;
    I2 = p.I2;
    g = p.g;

    z1 = q1 - q1_ref;
    z2 = dq1;

    v1 = -kp * z1 - kd * z2;

    

    d11 = m1*(lc1)^2 +m2*(l1^2+lc2^2+2*l1*lc2*cos(q2))+I1+I2;
    d12 = m2*(lc2^2+l1*lc2*cos(q2)) +I2;
    d21 = d12;
    d22 = m2*lc2^2+I2;

    
    h1 = -m2*l1*lc2*sin(q2)*dq2^2 - 2*m2*l1*lc2*sin(q2)*dq2*dq1;
    h2 = m2*l1*lc2*sin(q2)*dq1^2;

    
    phi1 = (m1*lc1+m2*l1)*g*cos(q1)+m2*lc2*g*cos(q1+q2);
    phi2 = m2*lc2*g*cos(q1+q2);

    ddq2 = -(d11 * v1 + h1 + phi1) / d12;

   
    
    if SIGMA == 1

        z1 = q1_ref - q1;
        z2 = dq1;
        v1 = kp * z1 - kd * z2;         % desired ddq1

        ddq2 = -(d11 * v1 + h1 + phi1) / d12;
        D1 =d21 -d22*d11/d12;
        H1 = h2 - h1*d22/d12;
        PHI1 = phi2 - phi1*d22/d12;
        u = D1*v1+PHI1+H1;
        %u = d21 * v1 + d22 * ddq2 + h2 + phi2;

    elseif SIGMA == 2
        D2 = d22-d21*d12/d11;
        H2 = h2 - h1 *d21/d11;
        PHI2 = phi2 - phi1*d21/d11;

      

        alpha = 0.7;
        eps = 0.2;
        sgn = atan(dq1 ); 
        % disp('SIgn')
        % disp(sgn);
        q2d   = alpha * sgn;

     

                v2   = -kp * (q2 - q2d) -kd *dq2;
        ddq1 = -(d12 * v2 + h1 + phi1) / d11;
        u =D2*v2+H2+PHI2;
        %u = d21 * ddq1 + d22 * v2 + h2 + phi2;

    
    end
   

    
   

end