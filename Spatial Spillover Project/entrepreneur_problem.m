function [policy_assets, value_function] = entrepreneur_problem(wages, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the entrepreneur's problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP PARAMETERS

% Value Function iteration parameters
max_iter = 2000; % Maximum number of iterations 
tol = 10^-2; % Not sensitive to lower values
Display = 0; % Display number of iterations

% Firm parameters
R = 1.04; % Intrest rate
beta = 0.85; % Discount rate
eta = params.eta; % Risk Aversion
n_villages = params.n_villages; % Number of villages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP GRIDS

% Set up grid for asset holdings.
asset_space = params.entrep_asset_space;
n_asset_states = length(asset_space);
asset_grid  = meshgrid(asset_space, asset_space);

% Set up the matricies for value function itteration
savind = zeros(n_asset_states, n_villages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UTILITY FUNCTION
% Pre-generate the period utility function

utility = zeros(n_asset_states, n_asset_states, n_villages);
profits = entrepreneur_profits(asset_space, wages, params);

% we have net_assets(i, j) = R * a_i  - a_j
net_assets = R.*asset_grid' - asset_grid;

% we have utility(i,j,k)  = utility of having assets i, and saving j, in village k
for zzz = 1:n_villages
    % Log utility 
    utility(:,:,zzz) = log(max(net_assets + repmat(profits(:,zzz), 1, n_asset_states), 10^-8));

    % CRRA utility 
    % consumption = net_assets + repmat(profits(:,zzz), 1, n_asset_states);
    % utility(:,:,zzz) = ((consumption.^(1-eta))-1)./(1-eta); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE VALUE FUNCTION
Vguess = squeeze(median(utility,2))./(1-beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITERATE ON VALUE FUNCTION

Vnext = Vguess;

Vdiff = 1;
iter = 0;

while iter <= max_iter && Vdiff > tol

    iter = iter + 1;
    Vlast = Vnext;

    % loop over villages
    for zzz = 1:n_villages

        Vchoice = bsxfun(@plus, utility(:,:,zzz), beta.*Vlast(:,zzz)');
        [Vnext(:,zzz), savind(:,zzz)] = max(Vchoice,[],2);
        p_assets(:,zzz) = asset_space(savind(:,zzz));

        % Adjust policy function for errors asssociated with log(0)
        for ii=(length(p_assets(:,zzz))-1):-1:1
            if p_assets(ii,zzz) > p_assets(ii+1,zzz)
                p_assets(ii,zzz) = 1;
            end
        end
    end

    Vdiff = norm(Vnext - Vlast, Inf);

    if Display >=1
        disp(['Entrepreneur Problem Iteration no. ' int2str(iter), ' max val fn diff is ' num2str(Vdiff)]);
    end
end

policy_assets = p_assets;
value_function = Vnext;

end
