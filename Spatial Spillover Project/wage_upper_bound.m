function [wage_ub] =  wage_upper_bound(params)

%% Solve For Wage Upper Bound
wage_ub = (1- params.alpha).*(params.z.*params.alpha./params.rho).^(params.alpha./(1- params.alpha));

end
