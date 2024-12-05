function [policy_assets, value_fun_out] = worker_problem(wages, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the households problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Parameters
% Value Function parameters
n_iterations = 2000; % Maximum number of iterations 
tol = 10^-2; % Not sensitive to lower values
% Household parameters
R = params.R; % Intrest rate
beta = params.beta; % Discount rate
kappa = params.kappa; % Migration cost matrix
a_bar = params.a_bar; % Borowing constraint
%Village Parameters 
n_villages  = params.n_villages; % Number of Villages 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up grid for asset holdings.

asset_space = params.worker_asset_space;
n_asset_states = length(asset_space);
asset_grid  = meshgrid(asset_space, asset_space);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Matrix of allowable moves for every asset level in every village
% An allowable move is one where a' - k_{ij} > - a_bar. That is, the budjet
% constraint must hold

trans_mat = false(n_villages, n_asset_states, n_villages);
asset_aft_move = zeros(n_villages, n_asset_states, n_villages);

for zzz = 1:n_villages
        for xxx = 1:n_villages
         trans_mat(zzz, : , xxx) = asset_space - kappa(zzz,xxx) >= - a_bar;
         asset_aft_move(zzz, : , xxx) = asset_space - kappa(zzz,xxx);
    end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the matricies for value function itteration
v_prime = zeros(n_asset_states, n_villages);
policy_assets = zeros(n_asset_states, n_villages); 

% We have net_assets(i, j) = R * a_i  - a_j
net_assets = R.*asset_grid' - asset_grid;

% Pre-generate the period utility function and feasibility conditions
for zzz = 1:n_villages
    % Log-Utility 
    utility(:,:,zzz) = log(max(net_assets + wages(zzz), 10^-8));
end

v_old = repmat(diag(median(utility,3)), 1, n_villages)./(1-beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commence value function itteration.

for iter = 1:n_iterations

    option_value = continuation_value(v_old, trans_mat, asset_aft_move, params);

    for zzz = 1:n_villages
    
        value_fun = bsxfun(@plus, utility(:,:,zzz), option_value(:,zzz)');
        [v_fun, ~ ] = max(value_fun,[],2);
        v_prime(:,zzz) = v_fun;

    end

    if norm(v_old - v_prime, Inf) < tol
        break
    end

    v_old = v_prime;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract policy functions

v_old = v_prime;
option_value = continuation_value(v_old, trans_mat, asset_aft_move, params);

for zzz = 1:n_villages

    value_fun = bsxfun(@plus, utility(:,:,zzz), option_value(:,zzz)');
    [v_fun, p_fun] = max(value_fun,[],2);

    policy_assets(:,zzz)= p_fun;
    value_fun_out(:,zzz) = v_fun;

end

end