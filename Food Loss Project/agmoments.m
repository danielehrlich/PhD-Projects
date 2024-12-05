function [moments] = agmoments(p,c,params)

    
    % % load parameters
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

    % aggregate moments
    Jm = (om/(om + phi)*pa-p)*x;
    theta = (Jm/kappa)^(1/alpha);
    ftheta = theta^(1-alpha);
    Vm = u(p*x);
    delta = d((M-c)/p_i);
    Vl = u(M);
    Vs = (u(c) + delta*Vl + ftheta* Vm)./(rho + ftheta + delta);
    fl = delta/(delta + ftheta);
    
    % save moments
    moments.p = p;
    moments.c = c;
    moments.i  = (M-c)/p_i;
    moments.theta = theta;
    moments.ftheta = ftheta;
    moments.delta = delta;
    moments.W = Vs;
    moments.fl = fl;

end