function muNew = evolveKinetic(mu,phaseL,phaseLOld,N,xi,pXi,TOL_ARNOLDI,M_ARNOLDI,dt,t,v,c,T,L) %

    xiL =  mInvL(m(xi,L),L);
    for i=1:1:2
      derPhaseXiL(i,:) = derivative_pr(xiL,phaseL,i,L);
      derPhaseXi(i,:) = trunk(derPhaseXiL(i,:));
    end
    phaseShNorm = -1j * derPhaseXi(1,:);
    phase = trunk(phaseL);
    phaseOld = trunk(phaseLOld);
    phaseDotNorm = -1j * (phase - phaseOld) / dt;
    termNorm = (phaseDotNorm - xi.*(Rd(t,v)/R(t,v)).*phaseShNorm);

    muNew = expmArnoldi(xi,mu,dt,v,TOL_ARNOLDI,M_ARNOLDI,c,pXi-calcA(t,T),pXi,derPhaseXi,termNorm,t,calcA(t,T),L);
end


