function [labor_demand, capital_demand, profits] = entrepreneur_input_choice(assets, wages, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
z = params.z; % TFP -  village level productivity
phi  = params.phi; % Leverage constraint
rho  = params.rho; % Rental rate of capital
alpha = params.alpha; % Alpha in production function
gamma = params.gamma; % Gamma in production function, Decreasing returns to scale if alapha + gamma < 1 
n_villages = params.n_villages; % Number of Villages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE INPUT CHOICE 

% Constant Returns to Scale
if alpha + gamma == 1 
	
	z_bar = rho/alpha * ((1-alpha)./wages).^((alpha - 1)./alpha);
	capital_demand = zeros(1, n_villages);
	max_k  = phi.*assets;
	capital_demand (z > z_bar) = max_k(z > z_bar);
	labor_demand = ((1-alpha)./wages).^(1/alpha).*z.*capital_demand;

% Decreasing Returns to Scale
elseif alpha + gamma < 1 

	capital_constrained = phi.*assets; 
	capital_unconstrained  = ((rho^(1-gamma)).*(alpha^(gamma-1))./z.*(gamma^(-gamma)).*(wages.^gamma)).^(1/(alpha + gamma -1));
	capital_demand = min(capital_constrained, capital_unconstrained);
	labor_demand = (wages./(gamma.*z.*(capital_demand).^(alpha))).^(1/(gamma-1));

% Increasing Returns to Scale
else
	disp('Error: Increasing Returns to Scale.');
end

% Calculate Profits
profits = z.*(capital_demand.^alpha).*(labor_demand.^gamma) - wages.*labor_demand - rho.*capital_demand;

% Check if profits > 0
labor_demand(profits < 0) = 0;
capital_demand(profits < 0) = 0;
profits(profits < 0) = 0;

end
