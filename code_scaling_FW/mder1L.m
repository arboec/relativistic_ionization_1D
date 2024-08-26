function mder1L = mder1L(xi,L)

    %mx = 110;%max(trunk(xi));
    mx = L;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    mder1L = a ./ (A*sqrt(1+a*a*xi.*xi));

end
