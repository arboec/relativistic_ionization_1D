function p = initializePGrid_pr(N,L)

  h = L/N;
  k = 0:N-1;
  p = k/(h*N);
  p = p - 0.5/h;
  p = p * 2*pi;
  p = p * 0.5; %because L -- is 1/2 from entire L

end
