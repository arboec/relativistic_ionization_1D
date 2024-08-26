function muNew = evolvePotential(t,v,dt,xi,muOld,LAMBDA,N) %lambda -- split operator parameter

        V = calcV(t,v,xi);
        VAbsorb = calcVAbsorb(xi,N);
        muNew = exp(-1j*dt*LAMBDA*(V+VAbsorb)).*muOld;
        %muNew = exp(-dt*LAMBDA*(V+0*VAbsorb)).*muOld; % for ground state

end
