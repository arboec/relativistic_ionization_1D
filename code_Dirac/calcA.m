function A = calcA(t,T)
        w = 2.0;
        T = 2*pi/w;
        %A = 0.3*(175/w) * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2); % correct to get Kylstra results
        %A =  2.0 * 0.3*(175/w) * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2); % correct to get Kylstra results
        %A = 3.8095 * 0.3*(175/w) * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2); % correct to get Kylstra results
        %this is what I showed before %A =  0.3*(175/w) * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2); % correct to get Kylstra results
        if (t > 2*pi) 
            A = 0;
        end

        %florescu, july
        w = 2;
        T = 2*pi/w;
        %A = 100 * sin(w*t) .* exp(-(1.1774*(t-2*T)/T).^2); %true florescu
        %modefied to make peaks sharper:
        %A = 100 * sin(w*t*0.8) .* exp(-(1.1774*(t-2*T)/T).^2) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2);
        A = 50 * sin(w*t) .* sin(2*t/(4*T)*pi/2).*sin(2*t/(4*T)*pi/2) .* (heaviside(2*T-t) + heaviside(t-2*T) .* exp(-(4*1.1774*(t-2*T)/T).^2));

%         %true florescue
%         w = 1;
%         T = 2*pi/w;
%         A = 50 * sin(w*t) .* exp(-(1.1774*(t-2*T)/T).^2); %true florescu
end
