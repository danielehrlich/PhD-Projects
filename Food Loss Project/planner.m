function [wmax, pmax, lmax] = planner(params)

    % load parameters
    rho     = params(1);
    gamma   = params(2);
    pa      = params(3);
    kappa   = params(4);
    x       = params(5);
    alpha   = params(6);
    d0      = params(7);
    z       = params(8);
    M       = params(9);
    abar    = params(10);
    beta    = params(11);
    p_i     = params(12);      

    % initialize utility and storage functions
    % eta =0.9;
    % u = @(x)((x.^(1-eta)-1)./(1-eta));
    u = @(x)(log(x-abar));
    u_prime = @(x)(1/(x-abar));
    d = @(x)(d0)*(1 + x).^(-beta);
    d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));

    % initialize search grid
    p_array = linspace(1e-4, pa, 100);
    l_array = linspace(1e-4, M, 20);
    [pgrid,lgrid] = meshgrid(p_array, l_array);

    % calculate welfare over grid
    delta = d((M-lgrid)/p_i);
    theta = [(pa-pgrid)*x/(rho + gamma)/kappa].^(1/alpha);
    ftheta = theta.^(1-alpha);
    s = gamma./(gamma +  ftheta + delta);
    v = 1/gamma .* ftheta .* s;
    y = 1/gamma .* delta .* s;
    h = (u(z + lgrid) + delta .*(u(z)/(rho + gamma)) + ftheta.*(u(z + pgrid*x)/(rho + gamma)))./(rho + (rho/(rho + gamma)) .* (delta + ftheta));
    w = 1/(rho + gamma) .* (u(z) + gamma*h);
    m = 1/(rho + gamma) .* (u(z + pgrid.*x)+ gamma*h);
    W = s.*h+ v.*m + y.*w;

    % find max
    wmax = max(W,[],"all");
    [lidx, pidx] = find(W == wmax);
    pmax = p_array(pidx);
    lmax = l_array(lidx);

end