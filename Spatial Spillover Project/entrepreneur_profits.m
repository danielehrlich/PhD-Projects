function [profits] = entrepreneur_profits(assets, wages, params) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS
z = params.z; % TFP -  village level productivity
phi  = params.phi; % Leverage constraint
rho  = params.rho; % Rental rate of capital
alpha = params.alpha; % Alpha in production function
gamma = params.gamma; % Gamma in production function, Decreasing returns to scale if alapha + gamma < 1 
n_villages = params.n_villages; % Number of Villages
entrep_na = params.entrep_na; % Asset spiace size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP GRID
profits = zeros(entrep_na, n_villages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE LABOR AND CAPITAL DEMAND
% Constant Returns to Scale
if alpha + gamma ==1
    
    for ia = 1:entrep_na
        profits(ia, :) = max(z.*alpha.*((1-alpha)./wages).^((1-alpha)/alpha)-rho, 0).*phi.*assets(ia);
    end

% Decreasing Returns to Scale
elseif alpha + gamma < 1 

	% Unconstrained Capital demand does not depend on assets and can be calculated outside of the loop
	capital_unconstrained  = ((rho^(1-gamma)).*(alpha^(gamma-1))./z.*(gamma^(-gamma)).*(wages.^gamma)).^(1/(alpha + gamma -1));

    for ia = 1:entrep_na
    	% Capital Demand given wages 
    	capital_constrained = phi.*assets(ia); 
		capital_demand = min(capital_constrained, capital_unconstrained);

		% Labor Demand given wages
		labor_demand = (wages./(gamma.*z.*(capital_demand).^(alpha))).^(1/(gamma-1));

		% Profits as a funciton of labor and capital demand 
        profits(ia, :) = z.*(capital_demand.^alpha).*(labor_demand.^gamma) - wages.*labor_demand - rho.*capital_demand;
	end

% Increasing Returns to Scale
else
	disp('Error: Increasing Returns to Scale.');
end

% If profits are less than 0, set profits to 0
profits(profits < 0) = 0;

end
