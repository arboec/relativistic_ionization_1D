function checkLinL = calcPhaseL(c,v,t,pExact,exactSolP,mult,xi,L)

  der1 = real(derivative_uniform_pr(pExact,exactSolP,1));

  check = exactSolP - pExact.*der1;
  xiCheck = trunk(der1);
  check = trunk(check);

  xi = m(xi,L);
  xiL = mInvL(xi,L);
  yL = xiL*R(t,v);
  RR = (yL.^16+c^16*t^16+1e-16).^(1/16); %1e-16 -- just a crutch against deivision by zero at t = 0
  x0L = c*t * mult*yL./RR;
  checkLinL = interp1(xiCheck,check,x0L,"spline","extrap"); %WARNING! SPLINE EXTRAPOLATION!

end
