function model_moments = agmoments(delta, lambda, s, g, Omega, params);
    
    % set parameters
    global qmax qlen r S qgrid_param;
    alpha = params(1);
    beta = params(2);
    gamma = params(3);
    phi = [1, params(4)];
    mbar = [params(5), 1-params(5)];
    theta = params(6);
    sigma = params(7);
    theta_hat = 1 - ((theta-1)/theta)^theta;
    phi_hat = [1 (phi(2)/phi(1))^theta; (phi(1)/phi(2))^theta 1];
    Profit = (1/(sigma-1))*(((sigma-1)/sigma)^sigma).*phi.^(1-sigma);
    qhat = linspace(0,1,qlen);
    qhat = qhat.^(1/qgrid_param);
    qhat = qmax*qhat;

    % probabilities 
    nu = 1 - sum(repmat(s.*delta,2,1).*alpha^theta.*(phi_hat.'),2)';  % [l, h]
    vfun = Profit./(1 - (beta/(beta+1))*(lambda.^(beta+1) + delta.^(beta+1)) - ((1-gamma)/(1+r)/((1+g)^(sigma-1)).*nu));
    psi = repmat(flip(s.*alpha^theta),2,1).*(phi_hat);
            
    % r&d share
    denom = (beta + 1).*(1 - (beta/(beta+1))*(lambda.^(beta+1) + delta.^(beta+1)) - ((1-gamma)/(1+r)/((1+g)^(sigma-1)).*nu));
    denom(denom<0) = 1e-6;
    rdshare = ((lambda.^(beta + 1) + delta.^(beta + 1))*(1/(beta + 1)))./(denom);
    rdshare(rdshare>1) = 1;
    agrdshare = rdshare*s';

    % probability of gaining/dropping a sector
    a = delta.*sum(psi,2)'; % probability of aquiring a product
    k = nu; % probability of keeping a product
    aga = a*s';
    agk = k*s';
    agc = 1 - agk; % creative destruction

    % employmnet growth rate
    empgi = (lambda.*(1+g)^(1-sigma)*(theta/(theta + 1 - sigma)) + (1-lambda).*((1+g)^(1-sigma))) - 1;
    empg = empgi*s';

    % relative scale 
    Opdf = s(1)*Omega(1,:) + s(2)*Omega(2,:);
    Ocdf = cumsum(Opdf);
    [~,scale01] = min(abs(Ocdf-0.999));
    rel_scale01 = (Opdf(1,scale01:end)*(qhat(scale01:end).^(sigma-1))'/sum(Opdf(1,scale01:end))/(Opdf*(qhat.^(sigma-1))'));

    % scope
    [xil, xih, M, E] = Xisolver(gamma, a, k, mbar);
    n_array  = 1:1:length(xil);
    Xipdf = s(1)*xil + s(2)*xih;
    Xicdf = cumsum(Xipdf);
    avscope = Xipdf*n_array';
    [~,scope01] = min(abs(Xicdf-0.999));
    rel_scope01 = (Xipdf(1,scope01:end)*n_array(scope01:end)'/sum(Xipdf(1,scope01:end))/avscope);
    [~,scope001] = min(abs(Xicdf-0.9999));
    rel_scope001 = (Xipdf(1,scope001:end)*n_array(scope001:end)'/sum(Xipdf(1,scope001:end))/avscope);
    
    
    % Exit rate
    agexit = gamma + E*M'/sum(M);
    
    % Model moments
    model_moments.g = g;
    model_moments.avscope = avscope;
    model_moments.rel_scope01 = rel_scope01;
    model_moments.rel_scope001 = rel_scope001;
    model_moments.rel_scale01 = rel_scale01;
    model_moments.a = a;
    model_moments.k = k;
    model_moments.aga = aga;
    model_moments.agk = agk;
    model_moments.agc = agc;
    model_moments.agexit = agexit;
    model_moments.empgi = empgi;
    model_moments.empg = empg;
    model_moments.delta = delta;
    model_moments.lambda = lambda;
    model_moments.g = g;
    model_moments.Omega = Omega;
    model_moments.s = s;

end

%%% OTHER POTENTIAL MOMENTS
% % growth rates 
% grl = zeros(1,nmax);
% grh = zeros(1,nmax);
% for qq = 1:nmax
%     suml = 0;
%     sumh = 0;
%     for xx = 1:nmax
%         suml = suml + n_array(xx) * tl(qq,xx);
%         sumh = sumh + n_array(xx) * th(qq,xx);
%     end
%     grl(qq) = suml/n_array(qq);
%     grh(qq) = sumh/n_array(qq);
% end
