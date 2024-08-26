function m = m(xi,L)

    %mx = 110;%max(xi);
    mx = L;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    m = asinh(xi*a)/A;

end
