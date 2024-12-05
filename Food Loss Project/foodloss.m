function [flmin, Wmin, pmin, lmin] = foodloss(params)

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
    eta =0.9;
    u = @(x)((x.^(1-eta)-1)./(1-eta));
    u_prime = @(x)(1/(x-abar));
    d = @(x)(d0)*(1 + x).^(-beta);
    d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));

    % initialize search grid
    p_array = linspace(1e-4, pa, 100);
    l_array = linspace(1e-4, M/p_i, 20);
    [pgrid,lgrid] = meshgrid(p_array, l_array);

    % calculate food loss over grid
    delta = d((M-lgrid)./p_i);
    theta = [(pa-pgrid)*x/(rho + gamma)/kappa].^(1/alpha);
    ftheta = theta.^(1-alpha);
    fl = delta./(delta + ftheta);
    output  = fl .* pgrid + theta * kappa + (M - lgrid)/p_i;

    % find min
    [cmin,I] = min(output,[],"all","linear");
    [lidx, pidx] = ind2sub(size(output),I);
    pmin = p_array(pidx);
    lmin = l_array(lidx);

    % calculate welfare
    delta = d((M-lmin)/p_i);
    theta = [(pa-pmin)*x/(rho + gamma)/kappa].^(1/alpha);
    ftheta = theta.^(1-alpha);
    flmin = delta/(delta + ftheta);
    s = gamma./(gamma +  ftheta + delta);
    v = 1/gamma .* ftheta .* s;
    y = 1/gamma .* delta .* s;
    h = (u(z + lmin) + delta .*(u(z)/(rho + gamma)) + ftheta.*(u(z + pmin*x)/(rho + gamma)))./(rho + (rho/(rho + gamma)) .* (delta + ftheta));
    w = 1/(rho + gamma) .* (u(z) + gamma*h);
    m = 1/(rho + gamma) .* (u(z + pmin.*x)+ gamma*w);
    Wmin = s.*h+ v.*m + y.*w;

end
