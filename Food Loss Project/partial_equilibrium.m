function [c, fl] = partial_equilibrium(params, p)

    % % load parameters
    % rho     = params(1);
    % kappa   = params(2);
    % x       = params(3);
    % alpha   = params(4);
    % d0      = params(5);
    % M       = params(6);
    % abar    = params(7);
    % beta    = params(8);
    % p_i     = params(9);    
    % eta     = params(10);
    % pa      = params(11);
    % % om      = params(12);
    % % phi     = 1/(1-om);

    % % initialize utility and storage functions
    % if eta == 1
    %     u = @(x)(log(x-abar));
    % else
    %     u = @(x)(((x-abar).^(1-eta)-1)./(1-eta));
    % end
    % d = @(x)(d0)*(1 + x).^(-beta);
    % d_prime = @(x)(-beta * d0*(1 + x).^(-beta-1));
    % % Omega = @(x)(min((omega * pa./x).^phi, 1));
    % Omega = @(x)(1-x/pa);

    % % initialize arrays
    % p_array = linspace(abar, pa, 100);
    % c_array = linspace(abar, M, 100);

    % % itterate on equilibrium until convergence
    % Vs = - 1;
    % iter = 0;
    % error = 1;
    % Vm_array = u(p_array.*x);
    % Jm_array = Omega(p_array).*1/2.*(pa - p_array)*x;
    % theta_array = (Jm_array./kappa).^(1/alpha);
    % ftheta_array = theta_array.^(1-alpha);
    % while ((iter < 300) && (error > 1e-15))
    %     Vhat = Vm_array - Vs;
    %     Vhat(Vhat < 0) = 0; 
    %     [~, pidx]  = max(ftheta_array.*Omega(p_array).*Vhat);
    %     p = p_array(pidx);
    %     % if Omega(p) < 1
    %     %     Jm = Omega(p).*x.*p./(phi-1);
    %     % else
    %     %     Jm = (pa * phi * om/(phi-1)-p)*x;
    %     % end
    %     Jm =  Omega(p)*1/2*(pa - p)*x;
    %     theta = (Jm/kappa)^(1/alpha);
    %     ftheta = theta^(1-alpha);
    %     Vm = u(p*x);
    %     Vs_array = (u(c_array) + ftheta*Omega(p) * Vm)./(rho + ftheta*Omega(p) + d((M-c_array)./p_i));
    %     [Vs_prime, cidx] = max(Vs_array);
    %     c = c_array(cidx);
    %     error = ((Vs-Vs_prime)/Vs)^2;
    %     iter = iter + 1;
    %     Vs = Vs_prime;
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
    Omega = @(x)(1 - betacdf(x, om, phi));

    % initialize arrays
    c_array = linspace(abar, M, 100);

    % itterate on equilibrium until convergence
    Vl = u(M);
    ex_value = (om / (om+ phi)) .* (1 - betacdf(p/pa, om+1, phi))/ (1 - betacdf(p/pa, om, phi));
    Jm =  Omega(p/pa)*(pa*ex_value - p)*x;
    theta = (Jm/kappa)^(1/alpha);
    ftheta = theta^(1-alpha);
    Vm = u(p*x);
    Vs_array = (u(c_array) +  d((M-c_array)./p_i).*Vl + ftheta*Omega(p/pa) * Vm)./(rho + ftheta*Omega(p/pa) + d((M-c_array)./p_i));
    [Vs_prime, cidx] = max(Vs_array);
    c = c_array(cidx);
    error = ((Vs-Vs_prime)/Vs)^2;
    iter = iter + 1;
    Vs = Vs_prime;
    fl = d((M-c)./p_i)/(d((M-c)./p_i) + ftheta*Omega(p/pa));

end