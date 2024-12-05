function [Omega] = Omegasolver(g, s, nu, psi, delta, lambda, params)

    % initialize params
    global qmax qlen r S qgrid_param ;
    alpha = params(1);
    beta = params(2);
    gamma = params(3);
    phi = [1, params(4)];
    mbar = [params(5), 1-params(5)];
    theta = params(6);
    sigma = params(7);
    theta_hat = 1 - ((theta-1)/theta)^theta;
    phi_hat = [1 (phi(2)/phi(1))^theta; (phi(1)/phi(2))^theta 1];
    Profit = (1/(sigma-1))*(((sigma-1)/sigma)^sigma).*phi.^(1-sigma);
    qhat = linspace(0,1,qlen);
    qhat = qhat.^(1/qgrid_param);
    qhat = qmax*qhat;

    lcdf = linspace(0, 1, qlen);
    hcdf = linspace(0, 1, qlen);

    for q = 1:(qlen)
        xidx(q) = max(find(qhat <=(1+g)*qhat(q)));
        qidx(q) = max(find(qhat <= (1+g)*qhat(q)/alpha));
    end

    % itterate on the distribution
    iter  = 0;
    error  = 1;
    while ((iter < 300) && (error > 1e-8))
        lpdf = [lcdf(1),diff(lcdf)];     
        hpdf = [hcdf(1),diff(hcdf)];
        for q=1:(qlen)
            lalpha = max(0, cumsum(lpdf.*(1 - (alpha./(((1+g)*qhat(q))./qhat)).^theta)));
            lsum = max(0, cumsum(lpdf.*(1 - (alpha./(((1+g)*qhat(q))./qhat)).^theta)));
            halpha = max(0, cumsum(hpdf.*(1 - (alpha./(((1+g)*qhat(q))./qhat)).^theta)));
            hsum = max(0, cumsum(hpdf.*(1 - (1./(((1+g)*qhat(q))./qhat)).^theta)));
            l1 = halpha(qidx(q));
            l2 = lalpha(xidx(q));
            l3 = lsum(xidx(q));
            h1 = lalpha(qidx(q));
            h2 = halpha(xidx(q));
            h3 = hsum(xidx(q));
            lcdf_prime(q) = 1/s(1)*((1-gamma)*(1-lambda(1))*nu(1)*s(1)*lcdf(xidx(q)) + ...
            (1-gamma)*s(1)*delta(1)*psi(1,2)*s(2)*l1 + (1-gamma)*s(1)*delta(1)*psi(1,1)*s(1)*l2 + ...
            (1-gamma)*lambda(1)*nu(1)*s(1)*l3+ gamma*mbar(1)*(s(1)*lcdf(q) + s(2)*hcdf(q)));
            hcdf_prime(q) = 1/s(2)*((1-gamma)*(1-lambda(2))*nu(2)*s(2)*hcdf(xidx(q)) + ...
            (1-gamma)*s(2)*delta(2)*psi(2,1)*s(1)*h1 + (1-gamma)*s(2)*delta(2)*psi(2,2)*s(2)*h2+ ...
            (1-gamma)*lambda(2)*nu(2)*s(2)*h3 + gamma*mbar(2)*(s(1)*lcdf(q)) + s(2)*hcdf(q));
        end 
        % lcdf_prime(qlen) = 1;
        % hcdf_prime(qlen) = 1;
        lcdf_prime = lcdf_prime/lcdf_prime(qlen);
        hcdf_prime = hcdf_prime/hcdf_prime(qlen);
        error = (lcdf_prime - lcdf)*(lcdf_prime - lcdf)';
        lcdf = lcdf_prime;
        hcdf = hcdf_prime;
        iter = iter + 1;
    end 

    lpdf = [lcdf(1),diff(lcdf)];
    hpdf = [hcdf(1),diff(hcdf)];
    Omega = [lpdf; hpdf];

end