function [moments] = agmoments_all(p,c,params)

    
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
    % om      = params(12);
    % phi     = 1/(1-om);
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

    % aggregate moments
    omega = Omega(p/pa);
    delta = d((M-c)/p_i);
    ex_value = (om / (om+ phi)) .* (1 - betacdf(p/pa, om+1, phi))/ (1 - betacdf(p/pa, om, phi));
    Jm =  Omega(p/pa)*(pa*ex_value - p)*x;
    theta = (Jm/kappa)^(1/alpha);
    ftheta = theta^(1-alpha);
    Vl = u(M);
    Vm = u(p*x);
    Vs = (u(c) + delta * Vl + ftheta * omega * Vm)./(rho + ftheta*omega + delta);
    fl = delta/(delta + ftheta*omega);
 
    % save moments
    moments.p = p;
    moments.c = c;
    moments.i  = (M-c)/p_i;
    moments.theta = theta;
    moments.ftheta = ftheta;
    moments.delta = delta;
    moments.W = Vs;
    moments.Vl = Vl;
    moments.fl = fl;
    moments.Omega = omega;

end