function muDer = derivative_pr(x,mu,order,L) 

    % derivatives are calculated via FFT with a "trick"
    % p * fP --> sin(2*p*h)/(2*h) * fP
    % such a substitute allows to calclate derivatives of non-periodic
    % functions at cost of lower accuracy (in fact, it is just finite differences)

    xNonUn = x;
    x = mL(x,L);
    h = (x(2)-x(1));
    N = size(x); N = N(2);
    p = initializePGrid_pr(N,0.5*h*N);

    muP = ifftshift(fft(fftshift(mu)));
    muDer1 = ifftshift(ifft(fftshift((1j*sin(2*p*h)/(2*h)) .* muP)));
    muDer = muDer1;
    if (order == 1)
        muDer = mder1L(xNonUn,L) .* muDer1;
    end
    if (order == 2)
        muDer = mder2L(xNonUn,L).*muDer1 + mder1L(xNonUn,L).*mder1L(xNonUn,L).*ifftshift(ifft(fftshift((1j*sin(2*p*h)/(2*h)).^2 .* muP)));
    end
end