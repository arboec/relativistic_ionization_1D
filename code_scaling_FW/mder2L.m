function mder2L = mder2L(xi,L)

    %mx = 110;%max(trunk(xi));
    mx = L;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    mder2L = -a*a*a*xi ./ (A*(1+a*a*xi.*xi).^(3/2));

end
