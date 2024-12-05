function [excess, results] = village_compute_soe(guess, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the small open economy version of Ehrlich and Townsend
% Given wages and parameters computer equilbirum for the economy
% Output difference in market clearing conditions and economy statistics

% 5 Major Steps:
% (1) Given wages, solve household problem
% (2) Given wages, solve entrepreneur's problem
% (3) Given policy functions from (1) and (2), computer invarient distribution
% (4) Check market clearing conditions
% (5) Compute economy statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wages = guess; % Wage guess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOUSEHOLD PROBLEM
% Solve the household's problem
[worker_policy_assets, worker_value_fun] = worker_problem(wages, params);

% Compute the invariant distribution for workers (across assets and locations) given
% the solution to the workers problem
worker_invariant_distribution = worker_invariant_dist(worker_value_fun, worker_policy_assets, params);

% Computer Worker Welfare and Consumption
[welfare_by_location, consumption_by_location] = worker_welfare(wages, params, worker_value_fun, worker_policy_assets, worker_invariant_distribution);

% Computer Labor Supply for every village
% Suming down a column gives the mass of workers in that location, i.e. labor supply
labor_supply = sum(worker_invariant_distribution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRM PROBLEM
[entrep_policy_assets, entrep_value_fun] = entrepreneur_problem(wages, params);

% Compute the invariant distribution for entrepreneurs (across assets and locations) given
% the solution to the entrepreneur's problem
entrep_invariant_distribution = entrepreneur_invariant_dist(entrep_policy_assets, params);

% Compute Labor Demand for every village, normalize by number of villages
[labor_demand, capital_demand, profit] = entrepreneur_input_choice(entrep_invariant_distribution, wages, params);
labor_demand = labor_demand./(params.pop_firm_ratio*params.n_villages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERROR IN WAGE GUESS
% Compute market clearing condition for labor
excess = labor_demand-labor_supply;

% Use Fishcer Burmeister Function to smooth multicomplementarity problem
upper_bound = (10.*ones(1,params.n_villages)).^10;
lower_bound = zeros(1,params.n_villages);
first_phi = fischer_burmeister((upper_bound) - (wages), excess);
excess = fischer_burmeister((wages) - (lower_bound), first_phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economy descriptive statistics

results.params = params;
results.wages = wages;
results.capital_demand = capital_demand;
results.labor_demand = labor_demand;
results.pop = labor_supply;
results.worker_policy_assets = worker_policy_assets;
results.entrep_policy_assets = entrep_policy_assets;
results.profit = profit;
results.worker_dist = worker_invariant_distribution;
results.entrep_dist = entrep_invariant_distribution;
results.worker_welfare = welfare_by_location;
results.consumption = consumption_by_location;

end
