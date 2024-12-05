function [invariant_distribution] = entrepreneur_invariant_dist(policy_assets, params)

%% SIMULATE FOR INVARIANT DISTRIBUTION)

% Parameters
Tsim = 500;

aspace = params.entrep_asset_space;
asim  = zeros(params.n_villages, Tsim);
[~,asim(:, 1)] = ismember(params.entrep_a0, aspace);
[~, pindex] = ismember(policy_assets, aspace);

%asset choice
for it = 1:(Tsim-1)
    for iv = 1:params.n_villages
        asim(iv, it+1) = pindex(asim(iv, it), iv);
    end
end

invariant_distribution = aspace(asim(:,Tsim));

end
