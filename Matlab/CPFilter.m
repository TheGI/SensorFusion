function [Quaternion, IntegralError] = CPFilter(q, a, w, m, dt, Kp, Ki, Ie)
  if(norm(a) == 0) 
    Quaternion = q;
    IntegralError = Ie;
    return
  end
  a = a / norm(a);

  if(norm(m) == 0) 
    Quaternion = q;
    IntegralError = Ie;
    return
  end
  m = m / norm(m);

  h = quaternProd(q, quaternProd([0 m], quaternConj(q)));
  b = [0 norm([h(2) h(3)]) 0 h(4)];

  x = [2*(q(2)*q(4) - q(1)*q(3)),
       2*(q(1)*q(2) + q(3)*q(4)),
       q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
  y = [2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3)),
       2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4)),
       2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2)]; 

  e = cross(a, x') + cross(m, y'); 
  if(Ki > 0)
    Ie = Ie + e * dt;
  else
    Ie = [0 0 0];
  end
  w = w + Kp * e + Ki * Ie;    
  qDot = 0.5 * quaternProd(q, [0, w(1), w(2), w(3)]);
  q = q + qDot * dt;
  Quaternion = q / norm(q);
  IntegralError = Ie;
endfunction