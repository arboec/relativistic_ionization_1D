function V = calcV(t,v,xi)

        %V = -exp(-R(t,v)*R(t,v)*xi.*xi/4);
        N = size(xi,2);
        EPS = 1e-3;
        V = (exp(-abs(R(t,v)*xi))-exp(-EPS*abs(R(t,v)*xi))) ./ abs(R(t,v)*xi);
        V(floor(0.5*N)) = (-1+EPS);
end
