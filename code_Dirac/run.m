% This is a MATLAB implementation of the code, written by Heiko Bauke in
% Python
% in "Computational Strong-Field Quantum Dynamics", edited by Dieter Bauer
% chapter 3 "Time-dependent relativistic wave equations: Numerics of the Dirac and the Klein-Gordon equation"
% listing 3, pages 106-107.


clear;

[N, L, SCALING, T, dt, c, GRIDUNI, FILENAME_RE, FILENAME_IM, TOL_ARNOLDI, M_ARNOLDI] = initializeSystemConstants();
[v,width,u] = initializeScaling(SCALING);

xIni = initializeGrid(N,L);
pIni = initializePGrid_pr(N,L);
Ninitial = N;
N = 1.5*N; L = 1.5*L;
x = initializeGrid(N,L);
p = initializePGrid_pr(N,L);

gInitial = initializeWaveFunction(FILENAME_RE,FILENAME_IM,Ninitial); gInitial = conj(gInitial');
g = interp1(xIni,gInitial,x,'spline');
g = g / trapz(x,abs(g));
gP = fftshift(fft(ifftshift(g)));

dt = 5e-5;
T = 2*pi;

x_init = -0.2;
p_mean = 106.4;
Delta_p=0.125*c;
p_0=sqrt(p.^2+c^2);
psi = zeros(2,N);
psi(1,:) = 1;
psi(2,:) = p./(c+p_0); 
psi = psi .* sqrt(0.5*(c+p_0)./p_0);
psiP = psi;
%gP = 1.0/(2*pi*Delta_p*Delta_p)^(0.25) * exp(-(p-p_mean).^2/4/Delta_p^2) .* exp(-1j*p*x_init);

psiP = psiP .* gP;
psiUpP = psiP(1,:); psiDownP = psiP(2,:);
psiUp = fftshift(ifft(ifftshift(psiUpP)));
psiDown = fftshift(ifft(ifftshift(psiDownP)));
psi(1,:) = psiUp; psi(2,:) = psiDown;

dPlus = sqrt(0.5*(c+p_0)./p_0);
Q = [diag(sparse(dPlus)) -diag(sparse(dPlus.*p./(c+p_0)))
     diag(sparse(dPlus.*p./(c+p_0))) diag(sparse(dPlus))];
QT = Q';
Umid = [diag(sparse(exp(-1j*dt*c*p_0))) 0*diag(sparse(exp(-1j*dt*c*p_0)))
        0*diag(sparse(exp(-1j*dt*c*p_0))) diag(sparse(exp(1j*dt*c*p_0)))];

rho2 = abs(psi(1,:)).^2 + abs(psi(2,:)).^2;
plot(x,abs(rho2));


counter = -1;
counterInternal = 1;
trajectoryX(500) = 0;
projecUp(500) = 0;
projecDown(500) = 0;
timePlot(500) = 0;

psiUp0 = psiUp;
t = 0;

tic
%=======
for t = 0:dt:9.425
V = calcV(x,N);
A = -calcA(t,T);

psiUp = psi(1,:); psiDown = psi(2,:);

psiUp = exp(-0.5*1j*V*dt).*psiUp;
psiDown = exp(-0.5*1j*V*dt).*psiDown;
psiUpInter = cos(-0.5*dt*c*abs(A))*psiUp + 1j*sign(A)*sin(-0.5*dt*c*abs(A))*psiDown;%sign(A)*
psiDown = 1j*sign(A)*sin(-0.5*dt*c*abs(A))*psiUp + cos(-0.5*dt*c*abs(A))*psiDown;%sign(A)*
psiUp = psiUpInter;
psi(1,:) = psiUp; psi(2,:) = psiDown;

psiUpP = fftshift(fft(fftshift(psiUp)));
psiDownP = fftshift(fft(fftshift(psiDown)));
psiP = Q*Umid*QT*conj([psiUpP, psiDownP]');
psiUpP = conj(psiP(1:N)');
psiDownP = conj(psiP(N+1:end)');
psiUp = fftshift(ifft(ifftshift(psiUpP)));
psiDown = fftshift(ifft(ifftshift(psiDownP)));

psi(1,:) = psiUp; psi(2,:) = psiDown;

psiUp = exp(-0.5*1j*V*dt).*psiUp;
psiDown = exp(-0.5*1j*V*dt).*psiDown;
psiUpInter = cos(-0.5*dt*c*abs(A))*psiUp + 1j*sign(A)*sin(-0.5*dt*c*abs(A))*psiDown;%sign(A)*
psiDown = 1j*sign(A)*sin(-0.5*dt*c*abs(A))*psiUp + cos(-0.5*dt*c*abs(A))*psiDown;%sign(A)*
psiUp = psiUpInter;

psi(1,:) = psiUp; psi(2,:) = psiDown;
counter = counter + 1;

if (mod(counter,5000) == 0)
    t
    rho2 = abs(psi(1,:)).^2 + abs(psi(2,:)).^2;
    trajectoryX(counterInternal) = trapz(x,x.*rho2)./trapz(x,rho2);
    
    psiUp = psi(1,:); psiDown = psi(2,:);
    psiUpP = fftshift(fft(fftshift(psiUp)));
    psiDownP = fftshift(fft(fftshift(psiDown)));
    rho2P = abs(psiUpP).^2 + abs(psiDownP).^2;
    projecUp(counterInternal) = trapz(p,abs(psiUpP).^2)./trapz(p,rho2P); %this is a share of upper component of wave function
    projecDown(counterInternal) = trapz(p,abs(psiDownP).^2)./trapz(p,rho2P); %this is a share of borrom upper component of wave function
    %these projections are sligthly different from projections computed in
    %the method code. In the method, the projections onto positive and
    %negative spinors are calculated, here we just compute share of upper
    %and bottom components The results are different by the factor of
    %"small component".
    
    timePlot(counterInternal) = t;
    counterInternal = counterInternal + 1;
        
    %------ transform back
    E = c*c*sqrt(1+p.^2/(c*c));
    U11 = (E+c*c)./ sqrt(2*E.*(E+c*c));
    U12 = (-p*c)./ sqrt(2*E.*(E+c*c));
    U21 = (p*c)./ sqrt(2*E.*(E+c*c));
    U22 = (E+c*c)./ sqrt(2*E.*(E+c*c));
    solP1b = U11 .* psiUpP + U21 .* psiDownP;
    solP2b = U12 .* psiUpP + U22 .* psiDownP;
    rPb = abs(solP1b).^2 + abs(solP2b).^2;
    normB(counterInternal)  = trapz(p,rPb);
    pS = p - A;
    E = c*c*sqrt(1+pS.^2/(c*c));
    U11 = (E+c*c)./ sqrt(2*E.*(E+c*c));
    U12 = (-pS*c)./ sqrt(2*E.*(E+c*c));
    U21 = (pS*c)./ sqrt(2*E.*(E+c*c));
    U22 = (E+c*c)./ sqrt(2*E.*(E+c*c));
    solP1b = U11 .* psiUpP + U21 .* psiDownP;
    solP2b = U12 .* psiUpP + U22 .* psiDownP;
    rPb = abs(solP1b).^2 + abs(solP2b).^2;
    normBS(counterInternal)  = trapz(p,rPb);
    %------    
end


end
%=======
timeOfcomp = toc

rho2 = abs(psi(1,:)).^2 + abs(psi(2,:)).^2;
rho2P = abs(psiUpP.^2 + abs(psiDownP).^2);
