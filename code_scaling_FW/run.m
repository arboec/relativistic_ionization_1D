N = 30720;
L = 400;
FILENAME_RE = sprintf('ground_state_kylstra_Re_L=400_shr_darwin_3N.dat'); % default ground state
FILENAME_IM = sprintf('ground_state_kylstra_Im_L=400_shr_darwin_3N.dat'); % default ground state

xiGS = initializeGrid(N,L); % mesh xi for ground state
pXiGS = initializePGrid_pr(N,L); % mesh of momentum pXi for ground state
muGS = initializeWaveFunction(FILENAME_RE,FILENAME_IM,N); muGS = conj(muGS');

[N, L, SCALING, T, dt, c, TOL_ARNOLDI, M_ARNOLDI] = initializeSystemConstants();
[v,width,u] = initializeScaling(SCALING);
xi = initializeGrid(N,L);
pXi = initializePGrid_pr(N,L);
xi = mInv(xi,L); %moving to physical coordinates instead of computtaional
mu = interp1(xiGS,muGS,xi,'spline'); %obtaining wave function for the actual mesh

pExact = initializeGrid(2*N,2*c); %neede to calculate phase
exactSolPA = 0;
phaseLOld = 0*expand(xi);

%=== initializing varibales...
t = 0; %starting simulation from time = 0;
counterVisual = 0;
counterProjection = 0;
counter = 1;
position(500) = 0;
Pminus(500) = 0;
Pplus(500) = 0;
time(500) = 0;
projection(500,2) = 0; 

multiplier = 1.05;%1.15
counterRefining = floor(timeOfNextRefining(multiplier,v)/dt)+1;
N0 = N;
%=== end of the initialization of varibales

tic
while (t <= 9.425)%4*T  8*T %3*pi 9.425

      if (counterVisual == counterRefining) % triggers mesh refinement
        t
        muP = fftshift(fft(ifftshift(mu)));
        M = 2*floor(1.05*N*0.5);
        muPNew = [ 0*[1:(M-N)/2] muP 0*[1:(M-N)/2]];
        pXi = initializePGrid_pr(M,L);
        muNew = (M/N) * fftshift(ifft(ifftshift(muPNew)));
        xiNew = initializeGrid(M,L);
        h = (xiNew(2)-xiNew(1));
        xiNew = xiNew+(M/N0-1)*h;

        xi = mInv(xiNew,L);
        mu = muNew;
        N = M;

        phaseLOld = calcPhaseL(c,v,t-dt,pExact,exactSolPA,u*(t-dt)/sqrt(1+u*u*(t-dt)*(t-dt)),xi,L);
        multiplier = multiplier * 1.05;
        tRefining = timeOfNextRefining(multiplier,v);
        counterRefining = floor(tRefining/dt)+1;
      end

    exactSolPA = calcExactSolutionA(exactSolPA,pExact,t,T,c,dt);
    phaseL = calcPhaseL(c,v,t,pExact,exactSolPA,u*t/sqrt(1+u*u*t*t),xi,L);

    muOld = mu;
    mu = evolvePotential(v,t,dt,xi,mu,0.5,N);
    mu = evolveKinetic(mu,phaseL,phaseLOld,N,xi,pXi,TOL_ARNOLDI,M_ARNOLDI,dt,t,v,c,T,L);
    mu = evolvePotential(v,t,dt,xi,mu,0.5,N);

    if (abs(mod(counterVisual,500)) == 0) % triggers computation of averages, such as <x> or share of positive/negative states
      t
      %spectrum_ravnomer(xi,phaseL,mu,t,v,pXi,L);
      rho = mu.*conj(mu);
      x = xi *R(t,v);
      counterProjection = counterProjection + 1;
      position(counterProjection) = trapz(x,x.*rho)/trapz(x,rho);
      projection(counterProjection,:) = calcProjections(xi,phaseL,mu,t,v,L,calcA(t,T),c);
      time(counterProjection) = t;
    end
    counterVisual = counterVisual + 1;

    phaseLOld = phaseL;
    t = t + dt; %increment while
end
timeOfComputation = toc

nameFile = ['data_dt=' num2str(dt) '_accuracy=' num2str(TOL_ARNOLDI) '_withProjections.mat'];
save(nameFile)
