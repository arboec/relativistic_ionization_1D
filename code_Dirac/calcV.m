function V = calcV(x,N)

        %c = 137;
        %V0 = 2*c^2+1e4;
        %w = 0.3/c;
        %V = 0.5*V0*(tanh(x/w)+1);
        %V = -exp(-x.*x*0.25);
       
        N = size(x,2);
        EPS = 1e-3;
        V = (exp(-abs(x))-exp(-EPS*abs(x))) ./ abs(x);
        V(floor(0.5*N)) = (-1+EPS);
        V = V + 1j*200*(heaviside(-(x-x(round(0.1*N)))).*(x-x(round(0.1*N))) - heaviside((x-x(N-round(0.1*N)))).*(x-x(N-round(0.1*N))));        
end
