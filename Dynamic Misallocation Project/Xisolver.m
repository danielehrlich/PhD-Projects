function [xil, xih, M, E] = Xisolver(gamma, a, k, mbar)

    global qmax qlen r S qgrid_param;

    %initialize arrays
    iterp = 0;
    nmax = 50;
    n_array  = 1:1:nmax;
    xil = ones(1,nmax)./nmax; % size distribution by type
    xih = ones(1,nmax)./nmax;
    Ml = S/2/(n_array*xil'); % number of firms by type
    Mh = S/2/(n_array*xih');

    % itterate on distribution
    while iterp < 30
        El = (1-gamma) * sum(xil .* ((1-a(1)).^n_array) .* ((1-k(1)).^n_array));
        Eh = (1-gamma) * sum(xih .* ((1-a(2)).^n_array) .* ((1-k(2)).^n_array));
        Ml_prime = Ml * ( 1- gamma - El) + mbar(1) * (gamma * S) ;
        Mh_prime = Mh * ( 1- gamma - Eh) + (1- mbar(1)) *  (gamma * S); 
        tl = zeros(nmax, n_array(nmax));
        th = zeros(nmax, n_array(nmax));
        for qq = 1:nmax
            for n = 1:n_array(nmax)
                q = n_array(qq);
                if q >= (n/2) & q <= n
                    suml = 0;
                    sumh = 0;
                    for p = (n-q):q
                        templ =  nchoosek(q,p) * nchoosek(q,n-p) * a(1)^p * (1-a(1))^(q-p) * k(1)^(n-p) * (1-k(1))^(q+p-n);
                        suml = suml + templ;
                        temph =  nchoosek(q,p) * nchoosek(q,n-p) * a(2)^p * (1-a(2))^(q-p) * k(2)^(n-p) * (1-k(2))^(q+p-n);
                        sumh = sumh + temph;
                    end
                    tl(qq,n) = (1-gamma) *suml;
                    th(qq,n) = (1-gamma) *sumh;
                end
                if q > n & q <= (2*n)
                    suml = 0;
                    sumh = 0;
                    for p = 0:n
                        templ = nchoosek(q,p) * nchoosek(q,n-p) * a(1)^p * (1-a(1))^(q-p) * k(1)^(n-p) * (1-k(1))^(q+p-n);
                        suml = suml + templ;
                        temph = nchoosek(q,p) * nchoosek(q,n-p) * a(2)^p * (1-a(2))^(q-p) * k(2)^(n-p) * (1-k(2))^(q+p-n);
                        sumh = sumh + temph;
                    end
                    tl(qq,n) = (1-gamma) *suml;
                    th(qq,n) = (1-gamma) *sumh;
                end        
                if q > (2*n)
                    tl(qq,n) = (1-gamma) * (1-a(1))^q * k(1)^n * (1-k(1))^(q-n);
                    th(qq,n) = (1-gamma) * (1-a(2))^q * k(2)^n * (1-k(2))^(q-n);
                end       
            end
        end
        % share of firms exiting
        tl0 = gamma + (1-gamma) .* ((1-a(1)).^n_array) .* ((1-k(1)).^n_array);
        th0 = gamma + (1-gamma) .* ((1-a(2)).^n_array) .* ((1-k(2)).^n_array); 
        % adjust upper bound to account for discretization
        tl(:, nmax) = tl(:, nmax) + (1 - sum([tl tl0'], 2));
        th(:, nmax) = th(:, nmax) + (1 - sum([th th0'], 2));
        tl( tl< 0 ) = 0;
        th( th< 0 ) = 0;
        
        % update disribution
        xil_prime = Ml / Ml_prime * sum(tl.*repmat(xil', [1,nmax]), 1);
        xih_prime = Mh / Mh_prime * sum(th.*repmat(xih', [1,nmax]), 1);
        xil_prime(1) = xil_prime(1) + mbar(1) * (gamma * S) / Ml_prime;
        xih_prime(1) = xih_prime(1) + (1-mbar(1)) * (gamma * S) / Mh_prime;
        xil_prime = xil_prime/sum(xil_prime);
        xih_prime = xih_prime/sum(xih_prime);
        xil = xil_prime;
        xih = xih_prime;
        Ml = Ml_prime;
        Mh = Mh_prime;
        iterp = iterp + 1;
        % if mod(iterp,100) == 0
        %     fprintf('At firm size iteration %d...\n',iterp);
        % end
    end
    E = [El, Eh];
    M = [Ml, Mh];

end