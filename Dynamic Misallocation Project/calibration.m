function [errorsq] = calibration(x0, data_moments)

    %% EQUILIBRIUM 
    global qmax qlen r S qgrid_param;
    params = x0;

    % change parameters if restrictions not met
    if params(1)*params(4) > 1
        params(4) = 1/params(2);
    end
    if (params(6) + 1 - params(7)) <= 0
        params(7) = params(6) + 0.9;
    end

    [delta, lambda, s, g, Omega, iter, error] = equilibrium(params);
    moments = agmoments(delta, lambda, s, g, Omega, params);

    % compare moments
    % model_moments = [moments.aga, moments.avscope,  moments.rel_scope01, moments.rel_scale01, moments.agexit, moments.g, moments.s(1)];
    model_moments = [moments.aga, moments.avscope,  moments.rel_scope01, moments.rel_scale01, moments.agexit, moments.g, moments.empg];
    error = (data_moments - model_moments)./data_moments;
    % increase the weight on relative scope and scale
    error(3) =  10*error(3);
    error(4) =  10*error(4);
    % error(6)  = 5* error(6);
    errorsq = error*error';

end