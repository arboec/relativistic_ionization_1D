function mInv = mInv(xi,L)

    %m = sqrt(1-xi.*xi);
    %mx = 110;%max(xi);
    mx = L;
    xi = xi/mx;
    %A = 4;
    A = 1.35;
    a = sinh(A)/mx;
    mInv = sinh(xi*A)/a;

end
