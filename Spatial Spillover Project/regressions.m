function [reg_data] = regression(old_ss, new_ss)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS

% Location of Villages
n = old_ss.params.n_villages;
t = old_ss.params.t_village; % Angle of each village along circle
R = 1; % Cirdle Radius
x0 = 0; % Center of the circle in the x direction.
y0 = 0; % Center of the circle in the y direction.
x = x0 + R.*cos(t); % x coordinates
y = y0 + R.*sin(t); % y coordinates

% Create buffer zone matrix
buffer_mat = zeros(old_ss.params.n_villages);
t_bar  = 4*pi/n;
for i = 1:n
    for j = 1:n
        if i == j
            buffer_mat(i,j) = 0;
        elseif abs(t(i)-t(j)) <= t_bar || abs(2*pi + t(i)-t(j)) <= t_bar || abs(t(i) - t(j) -2*pi) <= t_bar
            buffer_mat(i,j) = 1;
        end
    end
end


% Graphing Parameters
sz = 500; % Marker Size
fontSize = 15; % Font Size of Labels and Title

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables of interest

% Wages
wages1 = old_ss.wages; % Equilibirum Wages Before Shock
wages2 = new_ss.wages; % Equilibirum Wages After Shock
wage_diff = wages2-wages1; % Difference in Wages
log_wage_diff = log(wages2) - log(wages1); % Difference in log popualtion

% Population
pop1 = old_ss.pop; % Equilibirum population Before Shock
pop2 = new_ss.pop; % Equilibirum population After Shock
pop_diff = pop2-pop1; % Difference in population
log_pop_diff = log(pop2) - log(pop1); % Difference in log popualtion

% Leverage Constraint
lev1 = old_ss.params.phi; % Leverage constraint in first Equilibrium
lev2 = new_ss.params.phi; % Leverage constraint in second Equilibrium (after shock)
lev_diff = lev2-lev1; % Difference in Leverage Constraints
lev_av  = lev_diff*buffer_mat; % Average Leverage Difference of Neighbor

% Total Factor Productivity
tfp = old_ss.params.z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REGRESSION DATA

row_names = {'wage_diff', 'pop_diff', 'lev_diff', 'lev_av', 'tfp', 'log_wage_diff', 'log_pop_diff'}
reg_data = table(wage_diff', pop_diff', lev_diff', lev_av', tfp', log_wage_diff', log_pop_diff');
filename = 'reg_data.xlsx';
writetable(reg_data,filename,'Sheet',1,'Range','A1', 'WriteRowNames',true);

end
