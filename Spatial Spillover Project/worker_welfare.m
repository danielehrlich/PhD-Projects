function [welfare_by_location, consumption_by_location] = worker_welfare(wages, params, value_fun, policy_fun, invariant_distribution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE CONSUMPTION BY LOCATION

% Inputs:
% wages: Array of wage guesses for each village
% value_fun: Array of worker value functions for each location and each asset state
% policy_fun: Array of worker asset policy functions for each location and each asset state
% invariant_distribution: Worker invariant distribution matrix in steady state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
asset_grid = params.worker_asset_space;
n_asset_states = length(asset_grid);
n_villages = params.n_villages;
R = params.R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initailize arrays
assets_option = zeros(n_asset_states, n_villages);
consumption = zeros(n_asset_states, n_villages);
assets = zeros(n_asset_states, n_villages);
welfare_function = zeros(n_asset_states, n_villages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yyy = 1:n_asset_states
    for xxx = 1:n_villages

        assets(yyy,xxx) = asset_grid(policy_fun(yyy,xxx));
        consumption(yyy,xxx) = R.*asset_grid(yyy) + wages(xxx) - assets(yyy,xxx);
                        
    end
end

% Consumption
consumption_by_location = (sum(consumption.*invariant_distribution)./sum(invariant_distribution))';
agg_consumption = sum(sum(consumption.*invariant_distribution));

% Assets
assets_by_location = (sum(assets.*invariant_distribution)./sum(invariant_distribution))'; 
agg_assets = sum(sum(assets.*invariant_distribution));

% Welfare
welfare_by_location = (sum(value_fun.*invariant_distribution)./sum(invariant_distribution))';
agg_welfare = sum(sum(value_fun.*invariant_distribution));

end
