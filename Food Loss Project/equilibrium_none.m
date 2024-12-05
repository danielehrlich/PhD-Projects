function [c, p, iter, error] = equilibrium_none(params)

    % load parameters
    rho     = params(1);
    kappa   = params(2);
    x       = params(3);
    alpha   = params(4);
    d0      = params(5);
    M       = params(6);
    abar    = params(7);
    beta    = params(8);
    p_i     = params(9);    
    eta     = params(10);
    pa      = params(11);

    % initialize utility and storage functions
    if eta == 1
        u = @(x)(log(x-abar));
    else
        u = @(x)(((x-abar).^(1-eta)-1)./(1-eta));
    end
    d = @(x)(d0)*(1 + x).^(-beta);
    d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));

    % initialize arrays
    p_array = linspace(abar, pa, 100);
    c_array = linspace(abar, M, 100);

    % itterate on equilibrium until convergence
    Vs = 1;
    iter = 0;
    error = 1;
    Vl = u(M);
    Vm_array = u(p_array.*x);
    Jm_array = (pa-p_array)*x;
    theta_array = Jm_array./kappa;
    while ((iter < 300) && (error > 1e-15))
        Vhat = Vm_array - Vs;
        Vhat(Vhat < 0) = 0; 
        [~, pidx]  = max(((Vhat)).*(theta_array));
        p = p_array(pidx);
        Jm = (0.5*pa-p)*x;
        theta = (Jm/kappa);
        m = max(theta, 1);
        ftheta = m;
        Vm = u(p*x);
        Vs_array = (u(c_array) + Vl.* d((M-c_array)) + ftheta * Vm)./(rho + ftheta + d((M-c_array)./p_i));
        [Vs_prime, cidx] = max(Vs_array);
        c = c_array(cidx);
        error = ((Vs-Vs_prime)/Vs)^2;
        iter = iter + 1;
        Vs = Vs_prime;
    end
end