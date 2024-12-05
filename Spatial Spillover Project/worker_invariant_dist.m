function [invariant_distribution] = worker_invariant_dist(v_fun, policy_assets, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS

% Household parameters
R = params.R; % Intrest rate
beta = params.beta; % Discount rate
kappa = params.kappa; % Migration cost matrix
a_bar = params.a_bar; % Borowing constraint
nu = params.nu; % Elasticity of option value

%Village Parameters 
n_villages  = params.n_villages;
L = params.worker_dist; % initial distribution
asset_space = params.worker_asset_space;
n_asset_states = length(asset_space);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTRUCT MIGRATION SHARES
% want to construct a matrix with elements i,a,j, which is the probability that 
% worker in village i with assets a moves to village j

% Construct matrix of allowable moves between villages
trans_mat = false(n_villages, n_asset_states, n_villages);
asset_aft_move = zeros(n_villages, n_asset_states, n_villages);

for zzz = 1:n_villages
        for xxx = 1:n_villages
         trans_mat(zzz, : , xxx) = asset_space - kappa(zzz,xxx) >= - a_bar;
         asset_aft_move(zzz, : , xxx) = asset_space - kappa(zzz,xxx);
    end 
end 

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

% Apply migration share formula
v_fun_mat = exp((beta.*v_fun_mat)).^(1/nu);
v_fun_mat = fillmissing(v_fun_mat,'constant',0); 

% Construct migration share matrix. This the share of migrants 
% from i to j with asset state a'
migration_share = v_fun_mat./sum(v_fun_mat, 3);

% Change migration share matrix to share of migrants who start with asset state a
QM_mat =  zeros(n_villages, n_asset_states, n_villages);

for www = 1:n_villages
    for zzz=1:n_asset_states
        QM_mat(www,zzz,:) = migration_share(www, policy_assets(zzz, www), :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE ASSET TRANSITION MATRIX

% Probabiltiy that given that worker has assets a' at the end of the period
% conditional on living in village i, with assets a, and moving to village j
% Note that this probability does not take into account that someone in village i and assets a
% would never move to village j. We take care of this in the next step
QA_mat = zeros(n_villages, n_asset_states, n_villages, n_asset_states);

for zzz = 1:n_villages
    for www = 1:n_asset_states
        for xxx = 1:n_villages
            QA_mat(zzz,www,xxx,asset_idx(zzz,www,xxx)) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMBINE MIGRATION SHARES & ASSET TRANSITIONS
% Probabiltiy that worker living in village i, with assets a, moves
% to village j and has assets a' next period

asset_state_transition = zeros(n_asset_states, n_villages, n_asset_states, n_villages);


for www = 1:n_asset_states
    for zzz = 1:n_villages
        for xxx= 1:n_villages
           asset_state_transition(www, zzz, :, xxx) =  QA_mat(zzz, www, xxx, :).* QM_mat(zzz, www, xxx);
        end
    end
end

asset_state_transition = reshape(asset_state_transition,n_asset_states*n_villages, n_asset_states*n_villages);
asset_state_transition = sparse(asset_state_transition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITTERATE ON DISTRIBUTION 

L = reshape(L,1, n_asset_states*n_villages);

for zzz = 1:2000
    L_new = L*asset_state_transition;

    if norm(L_new-L) < 10^-10
        break
    end

    L = L_new;
end

invariant_distribution = reshape(L,n_asset_states, n_villages);


end
