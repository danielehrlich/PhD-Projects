function [c, p, iter, error] = equilibrium(params)

    % load parameters
    % rho     = params(1);
    % gamma   = params(2);
    % pa      = params(3);
    % kappa   = params(4);
    % x       = params(5);
    % alpha   = params(6);
    % d0      = params(7);
    % z       = params(8);
    % M       = params(9);
    % abar    = params(10);
    % beta    = params(11);
    % p_i     = params(12);    

    % % initialize utility and storage functions
    % % u = @(x)(log(x - abar));
    % % eta =0.9;
    % % u = @(x)((x.^(1-eta)-1)./(1-eta));
    % u = @(x)(log(x));
    % u_prime = @(x)(1/(x-abar));
    % d = @(x)(d0)*(1 + x).^(-beta);
    % d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));

    % % initialize arrays
    % p_array = linspace(1e-4, pa, 100);
    % l_array = linspace(1e-4, M, 20);

    % % itterate on equilibrium until convergence
    % l  = M-0.5;
    % ftheta = 27;
    % m_hat = 100;
    % error = 1;
    % iter = 0;
    % while ((iter < 1000) && (error > 1e-15))
    %     h = (rho + gamma)/(rho*(rho+ gamma + d((M-l)/p_i))) .* (d((M-l)/p_i)/(rho + gamma)*u(z) + u(z +l) +  ftheta*m_hat);
    %     m_array = 1/(rho + gamma) * (u(z + p_array.*x) - rho*h);
    %     m_array(m_array < 0) = 0; 
    %     j_array = 1/(rho + gamma)*(pa-p_array)*x;
    %     [~, pidx]  = max((m_array.^alpha).*(j_array.^(1-alpha)));
    %     p = p_array(pidx);
    %     j_hat =  1/(rho + gamma)*(pa-p)*x;
    %     theta = (j_hat/kappa)^(1/alpha);
    %     ftheta = theta^(1-alpha);
    %     m_hat =  1/(rho + gamma) * (u(z + p*x) - rho*h);
    %     h_array = (rho + gamma)./(rho*(rho+ gamma + d((M-l_array)./p_i))) .* (d((M-l_array)./p_i)./(rho + gamma).*u(z) + u(z + l_array) +  ftheta*m_hat);
    %     [h_prime, lidx] = max(h_array);
    %     l_prime = l_array(lidx);
    %     error = ((l - l_prime)/l)^2 + ((h-h_prime)/h)^2;
    %     iter = iter + 1;
    %     l = l_prime;
    % end

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
    om      = 1.3;
    phi     = 0.2;

    % initialize utility and storage functions
    if eta == 1
        u = @(x)(log(x-abar));
    else
        u = @(x)(((x-abar).^(1-eta)-1)./(1-eta));
    end
    d = @(x)(d0)*(1 + x).^(-beta);
    d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));
    % Omega = @(x)(min((omega * pa./x).^phi, 1));

    % initialize arrays
    p_array = linspace(abar+1, pa-1, 100);
    c_array = linspace(abar, M, 100);

    % itterate on equilibrium until convergence
    Vs = - 1;
    Vl = u(M);
    iter = 0;
    error = 1;
    Vm_array = u(p_array.*x);
    Jm_array = (om/(om + phi).*pa-p_array)*x;
    Jm_array(Jm_array<0) = 0;
    theta_array = (Jm_array./kappa).^(1/alpha);
    ftheta_array = theta_array.^(1-alpha);
    while ((iter < 300) && (error > 1e-15))
        Vhat = Vm_array - Vs;
        Vhat(Vhat < 0) = 0; 
        [~, pidx]  = max((Jm_array.^(1-alpha)).*(Vhat).^alpha);
        p = p_array(pidx);
        Jm = (om/(om + phi)*pa-p)*x;
        theta = (Jm/kappa)^(1/alpha);
        ftheta = theta^(1-alpha);
        Vm = u(p*x);
        Vs_array = (u(c_array) +  d((M-c_array)./p_i).*Vl + ftheta * Vm)./(rho + ftheta + d((M-c_array)./p_i));
        [Vs_prime, cidx] = max(Vs_array);
        c = c_array(cidx);
        error = ((Vs-Vs_prime)/Vs)^2;
        iter = iter + 1;
        Vs = Vs_prime;
    end
end