function A = calcA(t,T)
        %laser pulse from the article
        w = 2;
        T = 2*pi/w;
        A = 50 * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2) .* (heaviside(2*T-t) + heaviside(t-2*T) .* exp(-(4*1.1774*(t-2*T)/T).^2));
end
