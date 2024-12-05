function [errorsq] = calibration(x0, data_moments, params)

    alpha = x0(1);
    beta = x0(2);
    M = x0(3);
    % omega = x0(4);
    params(4)  = alpha;
    params(8)  = beta;
    params(6)  = M;
    % params(12) = omega;

    [c, p, iter, error] = equilibrium_all(params);
    [moments] = agmoments_all(p,c,params);
    fl = moments.fl;
    ftheta = moments.ftheta;
    delta_eq = moments.delta;
    % omega = moments.Omega;
    model_moments  = [ftheta, delta_eq, c];
    error  = (data_moments - model_moments)./data_moments;
    errorsq = error * error';

end