phase = trunk(phaseL-phaseL(1));

%   This file allows one to make a plots of:
%   1) plot(p, log_10(rho(p))) %  probability density in momentum space
%   2) plot(x, log_10(rho(x))) % probability density in coordinates
%   3) plot(time, position) % Trajectory
%   4) plot(time, projection(:,1-2)) % share of negative/positive states


% FFT interpolation:
% add point in P space
N = size(mu,2);
M = 2e5;
muP = fftshift(fft(ifftshift(mu)));
muPNew = [ 0*[1:(M-N)/2] muP 0*[1:(M-N)/2]];
pXiNew = initializePGrid_pr(M,L);
pNew = initializePGrid_pr(M,L*R(t,v));

% transfer back to coordinate space
muNew = M/N * fftshift(ifft(fftshift(muPNew)));
xiNew = initializeGrid(M,L);
xiNew = xiNew+(M/N0-1)*(xiNew(2)-xiNew(1));
xiNew = mInv(xiNew,L);

rhoNew = muNew .* conj(muNew);
plot(xiNew*R(t,v),log10(rhoNew/trapz(xiNew*R(t,v),rhoNew)),'.-')


%going to the uniform coordinate mesh
xiRavnomer = [-110:0.0015:110];
phaseRavnomer = interp1(xi,phase,xiRavnomer,'spline');
muRavnomer = interp1(xiNew,muNew,xiRavnomer,'spline');
%restoring the solution: multiply by cancelled phase
psiRavnomer = muRavnomer.*exp(-1j*phaseRavnomer);
psiRavnomerP = fftshift(fft(ifftshift(psiRavnomer)));
pNew = initializePGrid_pr(size(xiRavnomer,2),110*R(t,v));

rhoPRavnomer = psiRavnomerP.*conj(psiRavnomerP);

figure
plot(pNew,log10(rhoPRavnomer/trapz(pNew,rhoPRavnomer)),'.-')

%-----------------------------
%plot(time, -position)
%plot(time, projection(:,1), '-')