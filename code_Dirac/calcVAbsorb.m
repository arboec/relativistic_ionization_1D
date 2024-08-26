function VAbsorb = calcVAbsorb(xi,N)

            %N = size(xi); %i wish i had field .N for object xi...
            VAbsorb = 1j*200*(heaviside(-(xi-xi(round(0.1*N)))).*(xi-xi(round(0.1*N))) - heaviside((xi-xi(N-round(0.1*N)))).*(xi-xi(N-round(0.1*N))));

end
