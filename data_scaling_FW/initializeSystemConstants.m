function [N, L, SCALING, T, dt, c, TOL_ARNOLDI, M_ARNOLDI] = initializeSystemConstants()

    N = 7372;
    L = 110;
    SCALING = 1; %boolean, SCALING=1 -- "ON"; SCALING=0 -- "OFF"

    T = 2*pi; %time of calculations, a.u.
    dt = 0.0001; %step size along time line, a.u.

    c = 137; %spedd of light

    TOL_ARNOLDI = 1e-6;    %1e-8
    M_ARNOLDI = 45;         %100

end
