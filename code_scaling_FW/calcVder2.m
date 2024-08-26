function Vder2 = calcVder2(t,v,xi)

        N = size(xi,2);
        EPS = 1e-3;
        x = xi;
        Vder2 = exp(-abs(R(t,v)*x)).*(2+2*abs(R(t,v)*x)+R(t,v)*R(t,v)*x.*x)./abs(R(t,v)*R(t,v)*R(t,v)*x.*x.*x) - exp(-EPS*R(t,v)*abs(x)).*(2+2*EPS*R(t,v)*abs(x)+R(t,v)*R(t,v)*EPS*EPS*x.*x)./abs(R(t,v)*R(t,v)*R(t,v)*x.*x.*x);
        Vder2(floor(0.5*N)) = 1/3*(-1+EPS*EPS*EPS);
        Vder2 = Vder2 / (R(t,v)*R(t,v));
        
end
