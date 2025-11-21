function xL = expand(x)

       N = size(x); N = N(2);
       xL = [5*x(1)-4*x(2) 4*x(1)-3*x(2) 3*x(1)-2*x(2) 2*x(1)-x(2) x 2*x(N)-x(N-1) 3*x(N)-2*x(N-1) 4*x(N)-3*x(N-1) 5*x(N)-4*x(N-1)];

end
