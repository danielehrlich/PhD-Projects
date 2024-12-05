function [delta, lambda, s, g, Omega, iter, error] = equilibrium(params)

    % parameters
    global qmax qlen r S qgrid_param;
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

    % initialize values
    delta = [0.3 0.3];
    s = [0.5, 0.5];
    g = 0.15;
    Omega = [ones(1, qlen)/qlen; ones(1, qlen)/qlen];

    % parameter checks
    if (theta + 1 - sigma) <= 0
        disp("R1:  theta, sigma restriction not met.")
        sigma  = theta + 0.9;
    end
    if (alpha * phi(2)) > 1
        disp("R2: alpha, phi restriction not met")
        phi(2) = 1/alpha;
    end  

    %% SOLVE EQUATIONS
    iter = 0;
    error = 1;
    while ((iter < 100) && (error > 1e-4))

        % Update internal innovation rate
        nu= 1 - sum(repmat(s.*delta,2,1).*alpha^theta.*(phi_hat.'),2)'; 
        lambda = ((1-gamma)/(1+r)/((1+g)^(sigma-1))*(theta/(theta + 1 - sigma)-1).*nu).^(1/beta);  % [l, h]
        lambda = max(0, min(0.99, lambda));

        % Update type shares
        psi = [repmat(flip(s.*alpha^theta),2,1).*(phi_hat)];
        sl = ((1-gamma)*(delta(1)*psi(1,2)) + gamma * mbar(1))/((1-gamma)*(delta(2)*psi(2,1) + delta(1)*psi(1,2)) + gamma);
        s = [sl, 1-sl];  % [l, h]

        % calculate growth rate
        type_g(1) = 1 + lambda(1) * nu(1) * (theta/(theta + 1 -sigma)-1) + s(1)*delta(1)*alpha^theta*(alpha* theta/(theta+1-sigma)-1) + s(2)*delta(2)*phi_hat(2,1)*alpha^theta*((phi(1)/phi(2))*alpha*theta/(theta+1-sigma)-1);
        type_g(2) = 1 + lambda(2) * nu(2) * (theta/(theta + 1 -sigma)-1) + s(2)*delta(2)*alpha^theta*(alpha* theta/(theta+1-sigma)-1) + s(1)*delta(1)*phi_hat(1,2)*alpha^theta*((phi(2)/phi(1))*alpha*theta/(theta+1-sigma)-1);
        type_qpsi = s./(phi.^(sigma-1)).*(sum(Omega.*repmat(qhat, 2, 1), 2))';
        qpsi = sum(type_qpsi);
        g = ((type_g * type_qpsi')./qpsi)^(1/(sigma-1))-1;
        g = max(g,0);

        % Update product quality distributions
        Omega = Omegasolver(g, s, nu, psi, delta, lambda, params);
    
        % Update external innovation rates
        vfun = Profit./(1 - (beta/(beta+1))*(lambda.^(beta+1) + delta.^(beta+1)) - ((1-gamma)/((1+r)*((1+g)^(sigma-1))).*nu));
        vfun(vfun<0) = 1e6;
        avqual = [(Omega(1,:)*qhat').^(sigma-1) (Omega(2,:)*qhat').^(sigma-1)];
        vbar = kron(vfun',avqual); % [ll lh; hl hh]
        delta_prime = ((1-gamma)/(1+r)*((alpha^theta * theta)/(theta + 1 - sigma))./((1+g).^(sigma-1)).*sum(psi .* vbar, 2)).^(1/beta);
        delta_prime = max(min(delta_prime, 1), 0);

        % Update distributions
        error = ((delta_prime' - delta)./delta)*((delta_prime' - delta)./delta)';
        delta = delta_prime';
        iter = iter + 1;
    end

end

