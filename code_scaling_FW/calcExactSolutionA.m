function exactSolPA = calcExactSolutionA(exactSolPA,p,t,T,c,dt) %

    A0 = calcA(t,T);
    exactSolPA = exactSolPA + c*c*dt*sqrt(1+((p-A0)/c).^2);

end

