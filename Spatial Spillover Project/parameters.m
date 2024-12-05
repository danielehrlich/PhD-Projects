function [params] = parameters()

% Draw random Numbers
rng(2020);

% Household parameters
params.eta = 0.9; % Risk-Aversion
params.beta = 0.926; % Discount factor
params.R =  1.054; %1/params.beta; % Interest rate on assets
params.nu = 1/3; % Contintuation value elasticity

% Asset grid parameters
params.a_bar = -0.08;
params.worker_asset_space = clustergrid(params.a_bar, 0.50, 6, 20, 30, 0.7, 0.3);
params.worker_na = length(params.worker_asset_space);
params.entrep_asset_space = clustergrid(1, 2, 20, 20, 30, 0.7, 0.3);
params.entrep_na = length(params.entrep_asset_space);

% Village Parameters
params.n_villages = 10; % Number of village
params.beta_a = 1; % First paramter of Beta Distribution
params.beta_b = 9; % Second paramter of Beta Distribution
[params.t_village, params.village_dist] = village_locations_circle(params.n_villages, params.beta_a, params.beta_b);
params.kappa = 1*params.village_dist + 0*params.village_dist^2; %Migration cost matrix 
params.worker_dist = ones(params.worker_na, params.n_villages)/(params.worker_na *params.n_villages);
params.pop_firm_ratio = 40;

% Entrepreneur Parameters
params.entrep_beta = 0.9; % Discount factor
params.z = 0.96.*randn(1, params.n_villages) + 3.74; % TFP -  village level productivity
params.z(params.z <= 0) = 1; % tfp should be bounded at 0. For good measure set to 1.
params.phi  = 1.2.*ones(1, params.n_villages); % Leverage constraint
params.rho  = .04; % Rental rate of capital
params.entrep_a0 = params.entrep_asset_space(5).*ones(1, params.n_villages); % Initial asset endownment of entrepreneurs
params.alpha = 0.27;        
params.gamma = 0.16;

end
