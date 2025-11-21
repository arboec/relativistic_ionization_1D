function VAbsorb = calcVAbsorb(xi,N)

    VAbsorb = 1j*10*(heaviside(-(xi-xi(round(0.045*N)))).*(xi-xi(round(0.045*N))) - heaviside((xi-xi(N-round(0.045*N)))).*(xi-xi(N-round(0.045*N))));

end
