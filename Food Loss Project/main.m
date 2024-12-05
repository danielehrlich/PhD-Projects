
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho     = 0.03;
kappa   = 10;
x       = 1;
alpha   = 0.3;
d0      = 20;
M       = 15;
abar    = 1.1;
beta    = 0.35;
p_i     = 1;    
eta     = 1;
pa      = 770;
% omega   = 0.2;
% phi     = 10;
params = [rho, kappa, x, alpha, d0, M, abar, beta, p_i, eta, pa];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [alpha, beta, M, kappa];
data_moments = [22, 2, 10, 0.9];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1, 0.1, abar, 0];
ub = [1, 2, 100, 100];
fun = @(x)calibration_all(x, data_moments, params);
options = optimoptions('fmincon','Display','iter','MaxIterations', 50);
[elas,fval,exitflag,output] =  fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
params(4)  = elas(1);
params(8)  = elas(2);
params(6)  = elas(3);
params(2)  = elas(4);
save paramsmat params;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Gains from reducing frictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline with all frictions
load paramsmat;
[ceq, peq, iter, error] = equilibrium_all(params);
[moments_eq_all] = agmoments_all(peq,ceq,params);
c_equiv_eq_all = abar + exp(moments_eq_all.W);
fl_eq_all  = moments_eq_all.fl;
peq_all = peq;
ceq_all = ceq;

% without limitted commitment but search 
load paramsmat;
[ceq, peq, iter, error] = equilibrium(params);
[moments_eq_lc] = agmoments(peq,ceq,params);
c_equiv_eq_lc = abar + exp(moments_eq_lc.W);
fl_eq_lc  = moments_eq_lc.fl;
peq_lc = peq;
ceq_lc = ceq;

% without search but limited commitment 
load paramsmat;
params(4) = 0.1;
[ceq, peq, iter, error] = equilibrium_all(params);
[moments_eq_all1] = agmoments_all(peq,ceq,params);
c_equiv_eq_all1 = abar + exp(moments_eq_all1.W);
fl_eq_all1  = moments_eq_all1.fl;
peq_all1 = peq;
ceq_all1 = ceq;

% without limitted commitment or search 
load paramsmat;
params(4) = 0.1;
[ceq, peq, iter, error] = equilibrium(params);
[moments_eq_lc2] = agmoments(peq,ceq,params);
c_equiv_eq_lc2 = abar + exp(moments_eq_lc2.W);
fl_eq_lc1  = moments_eq_lc2.fl;
peq_lc1 = peq;
ceq_lc1 = ceq;

c_equiv = [1, c_equiv_eq_lc/c_equiv_eq_all, c_equiv_eq_all1 /c_equiv_eq_all, c_equiv_eq_lc2 /c_equiv_eq_all];
fl = [fl_eq_all, fl_eq_lc, fl_eq_all1, fl_eq_lc1];
p_eq = [peq_all, peq_lc, peq_all1, peq_lc1];
c_eq = [ceq_all, ceq_lc, ceq_all1, ceq_lc1];
table(c_equiv', fl', p_eq', c_eq')
