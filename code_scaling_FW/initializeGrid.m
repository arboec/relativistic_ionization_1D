function xi = initializeGrid(N,L)

    h = 2 * L / N;
    xi(N) = 0;
    for i=1:1:N
        xi(i) = -L + h*i;
    end
    
end