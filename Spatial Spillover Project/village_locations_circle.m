function [t, village_dist] = village_locations_circle(n, alpha, beta)

    % Create a random set of coordinates in a circle.
    % Create matrix of allowable location moves for workers

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PARAMETERS
    R = 1/pi; % Cirdle Radius
    x0 = 0; % Center of the circle in the x direction.
    y0 = 0; % Center of the circle in the y direction.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CREATE SET OF POINTS
    
    % Generate random numbers from beta distribution
    draws = betarnd(alpha,beta,n,1);
    
    % For a full circle, use 0 and 2*pi.
    angle1 = 0;
    angle2 = 2*pi;
    t = (angle2 - angle1) * draws + angle1;
    t = sort(t);
    
    % Calculate distance between every village
    village_dist = zeros(n);
   
    for i = 1:n
        for j = 1:n
            village_dist(i,j) = R * min(2*pi - mod(t(i) - t(j), 2*pi), mod(t(i) - t(j), 2*pi));
        end
    end

end
