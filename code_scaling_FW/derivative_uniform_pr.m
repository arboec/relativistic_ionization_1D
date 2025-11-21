function muDer = derivative_uniform_pr(x,mu,order) %

    % derivatives are calculated via FFT with a "trick"
    % p * fP --> sin(2*p*h)/(2*h) * fP
    % such a substitute allows to calclate derivatives of non-periodic
    % functions at cost of lower accuracy (in fact, it is just finite differences)

    h = (x(2)-x(1));
    N = size(x); N = N(2);
    p = initializePGrid_pr(N,0.5*h*N);

    muP = ifftshift(fft(fftshift(mu)));
    muDer = ifftshift(ifft(fftshift((1j*sin(2*p*h)/(2*h)).^order .* muP)));
    
end

