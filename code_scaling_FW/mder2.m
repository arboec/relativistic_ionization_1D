function mder2 = mder2(xi,L)

    %mx = 110;%max(xi);
    mx = L;
    A = 1.35;
    a = sinh(A)/mx;
    A = A/mx;
    mder2 = -a*a*a*xi ./ (A*(1+a*a*xi.*xi).^(3/2));

end
