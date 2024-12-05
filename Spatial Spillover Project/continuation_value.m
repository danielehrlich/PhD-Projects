function [option_value] = continuation_value(v_fun, trans_mat, asset_aft_move, params)

% Parameters
n_villages = params.n_villages;
asset_space = params.worker_asset_space;
n_asset_states = length(asset_space);
kappa = params.kappa;
nu = params.nu;
beta = params.beta;

% Initiate Arrays
option_value = zeros(n_asset_states, n_villages);

% Find Closest element in asset space to assets after migration
asset_aft_move = reshape(asset_aft_move, n_villages*n_asset_states*n_villages, 1);
edges = [-Inf, mean([asset_space(2:end); asset_space(1:end-1)]), +Inf];
asset_idx = discretize(asset_aft_move, edges);
asset_idx = reshape(asset_idx, n_villages, n_asset_states, n_villages);

% Create linear index of villages
village_idx = find(asset_idx);
[~, ~, village_idx] = ind2sub(size(asset_idx), village_idx);
village_idx = reshape(village_idx, n_villages, n_asset_states, n_villages);

% Create index of both villages and assets
idx  = sub2ind(size(v_fun), asset_idx, village_idx);

% Extract value function elements associated with index
v_fun_mat = v_fun(idx);
v_fun_mat = v_fun_mat.*trans_mat;
drop_idx = find(abs(v_fun_mat) == 0);
v_fun_mat(drop_idx) = NaN;

% Compute Option Value
for zzz = 1:n_villages
    for www = 1:n_asset_states
    	
    	summation = nansum(exp(beta.*v_fun_mat(zzz,www,:)).^(1/nu));
        option_value(www,zzz) = nu * log(summation);

    end
end

end