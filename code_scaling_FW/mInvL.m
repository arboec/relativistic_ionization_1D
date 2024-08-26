function mInvL = mInvL(xi,L)

    %mx = 110;%max(xi);
    mx = L;
    xi = expand(xi);
    xi = xi/mx;
    A = 1.35;
    a = sinh(A)/mx;
    mInvL = sinh(xi*A)/a;

end
