function [G,r] = ct_zero(x,y)
%CT_ZERO Compute the core transformation that introduces a zero in y, could be replace by planerot

  if (y == 0)
  c = 1; s = 0; r = x;
  else
    if (abs(x) >= abs(y))
      theta=conj(sign(x));
      t = y/x; r = sqrt(1 + abs(t)^2);
      c = theta/r;
      s = conj(t)*c;
      r = theta*x*r;
    else
      theta=conj(sign(y));
      t = x/y; r = sqrt(1 + abs(t)^2 );
      s = theta/r;
      c = conj(t)*s;
      r = theta*y*r;
    end
  end
  G = [c s; -conj(s) conj(c)];
end