function projections = calcProjections(xi,phaseL,mu,t,v,L,A,c)

phase = trunk(phaseL-phaseL(1));

N = size(mu,2);
M = 2e5;
muP = fftshift(fft(ifftshift(mu)));
muPNew = [ 0*[1:(M-N)/2] muP 0*[1:(M-N)/2]];

muNew = M/N * fftshift(ifft(fftshift(muPNew)));
xiNew = initializeGrid(M,L);
xiNew = xiNew+(M/N-1)*(xiNew(2)-xiNew(1));
xiNew = mInv(xiNew,L);



xiRavnomer = [-110:0.0015:110];
phaseRavnomer = interp1(xi,phase,xiRavnomer,'spline');
muRavnomer = interp1(xiNew,muNew,xiRavnomer,'spline');
psiRavnomer = muRavnomer.*exp(-1j*phaseRavnomer);
psiRavnomerP = fftshift(fft(ifftshift(psiRavnomer)));
pNew = initializePGrid_pr(size(xiRavnomer,2),110*R(t,v));


p = pNew - A;
E = c*c*sqrt(1+p.^2/(c*c));
U11 = (E+c*c)./ sqrt(2*E.*(E+c*c));
U12 = (-p*c)./ sqrt(2*E.*(E+c*c));
U21 = (p*c)./ sqrt(2*E.*(E+c*c));
U22 = (E+c*c)./ sqrt(2*E.*(E+c*c));

up = U11.*psiRavnomerP;
down = U21.*psiRavnomerP;

p = pNew;
normal = trapz(p,abs(up).^2 + abs(down).^2);
E = c*c*sqrt(1+p.^2/(c*c));

U11 = (E+c*c)./ sqrt(2*E.*(E+c*c));
U12 = (-p*c)./ sqrt(2*E.*(E+c*c));
U21 = (p*c)./ sqrt(2*E.*(E+c*c));
U22 = (E+c*c)./ sqrt(2*E.*(E+c*c));

solP1c = (U11 .* up + U21 .* down) / sqrt(normal);
solP2c = (U12 .* up + U22 .* down) / sqrt(normal);
projections(1) = trapz(p,abs(solP1c).^2);
projections(2) = trapz(p,abs(solP2c).^2);

end