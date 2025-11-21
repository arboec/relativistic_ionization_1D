function mder1 = mder1(xi,L)

    %mx = 110;%max(xi);
    mx = L;
    %A = 4;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    mder1 = a ./ (A*sqrt(1+a*a*xi.*xi));

end
