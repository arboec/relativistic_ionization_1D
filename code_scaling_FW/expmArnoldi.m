%   This is a MATLAB implementation of code, originally written 
%   by Mike A. Botchev in "A short guide to exponential Krylov subspace 
%   time integration for Maxwell s equations", Figure 1, page 11, ISSN
%   1874-4850, available from: www.math.utwente.nl/publications.
%   

function y = expmArnoldi(xi,mu,dt,v,toler,m,c,pXi,pXi0Field,derPhaseANum,termNorm,t,A0,L)
% y = expm_Arnoldi(A,v,t,toler,m)
% computes $y = \exp(-t A) v$
% input: A (n x n)-matrix, v n-vector, t>0 time interval,
% toler>0 tolerance, m maximal Krylov dimension

mu = conj(mu');

H = zeros(m+1,m);

beta = norm(mu);
V = repmat(mu/beta,1,m);
V(:,1) = mu/beta;
term = calculateCoefficients(R(t,v),A0,c,derPhaseANum(1,:),derPhaseANum(2,:)/(R(t,v)*R(t,v)));


for j=1:m
    w = applyOperator(t,v,xi,c,conj(V(:,j))',pXi0Field,term,L);
    w = conj(w');

    for i=1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end

    H(j+1,j) = norm(w);
    e1 = zeros(j,1); e1(1) = 1;
    ej = zeros(j,1); ej(j) = 1;
    s = [0.01, 1/3, 2/3, 1]*dt;

    for q=1:length(s)
        u = expm(-s(q)*H(1:j,1:j))*e1;
        beta_j(q) = -H(j+1,j)* (ej'*u);
    end

    resnorm = norm(beta_j,'inf');
    %fprintf('j = %d, resnorm = %.2e\n',j,resnorm);
    if resnorm<=toler
        break
        elseif j==m
        disp('warning: no convergence within m steps');
    end

    V(:,j+1) = w/H(j+1,j);
end

y = V(:,1:j)*(beta*u);
y = conj(y)';
y = exp(-dt*(1j*c*c * term(1,:) + termNorm)) .* y;
end
