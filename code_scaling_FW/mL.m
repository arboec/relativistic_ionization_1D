function mL = mL(xi,L)

    %mx = 110;%max(trunk(xi));
    mx = L;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    mL = asinh(xi*a)/A;

end
