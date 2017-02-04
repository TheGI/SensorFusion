function [Quaternion, BiasError] = GDFilter(q, g, w, m, dt, beta, zeta, w_b)
  if(norm(g) == 0) 
    Quaternion = q;
    BiasError = w_b;
    return
  end
  g = g / norm(g);

  
  %if(norm(m) == 0)
  %  Quaternion = q;
  %  BiasError = w_b;
  %  return
  %end
  %m = m / norm(m);
  %h = quaternProd(q, quaternProd([0 m], quaternConj(q)));
  %b = [0 norm([h(2) h(3)]) 0 h(4)];
  

  f = [2*(q(2)*q(4) - q(1)*q(3)) - g(1);
                2*(q(1)*q(2) + q(3)*q(4)) - g(2);
                2*(0.5 - q(2)^2 - q(3)^2) - g(3)];
                
                %2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3)) - m(1);
                %2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4)) - m(2);
                %2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2) - m(3)];
                
  J = [-2*q(3), 2*q(4), -2*q(1), 2*q(2);
       2*q(2), 2*q(1), 2*q(4), 2*q(3);
       0, -4*q(2), -4*q(3), 0];
       %-2*b(4)*q(3), 2*b(4)*q(4), -4*b(2)*q(3)-2*b(4)*q(1), -4*b(2)*q(4)+2*b(4)*q(2);
       %-2*b(2)*q(4)+2*b(4)*q(2), 2*b(2)*q(3)+2*b(4)*q(1),	2*b(2)*q(2)+2*b(4)*q(4), -2*b(2)*q(1)+2*b(4)*q(3);
       %2*b(2)*q(3), 2*b(2)*q(4)-4*b(4)*q(2), 2*b(2)*q(1)-4*b(4)*q(3), 2*b(2)*q(2)];
       
  % Quaternion Error Estimation
  dq = J'*f;
  dq = (dq / norm(dq))';
  
    
  % Gyro Bias Correction
  w_e = 2 * quaternProd(quaternConj(q), dq);
  w_b = w_b + w_e(2:4)*dt*zeta;
  w = w - w_b;

  % Update new Quaternion
  dq_w = 0.5 * quaternProd(q, [0 w]);
  q = q + (dq_w - beta*dq)*dt;
  q = q / norm(q);
  
  Quaternion = q;
  BiasError = w_b;
end
