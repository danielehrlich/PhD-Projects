clear all;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set global variables
rng default;
global qmax qlen r S qgrid_param;
qlen = 50;
qmax = 2;
r = 0.3;
S = 10000;
qgrid_param = 0.3; %1 for linear, 0 for L-shaped
% sigma = 3.5;

% agregate moments  
g = .15;
aga = 0.3;
agk = 0.6;
empg = 0.03;
avscope = 1.1;
rel_scale01 = 53;
rel_scope01 = 7.7;
s = 0.01;
agexit = 0.6;
data_moments = [aga, avscope, rel_scope01, rel_scale01, agexit, g, empg];

% % parameter guess
alpha = 0.81;
beta = 15.9;
gamma = 0.1;
phih =  1.02;
mbarl = 0.1;
theta =  8;
sigma = 4;
x0 = [alpha, beta, gamma, phih, mbarl, theta, sigma];
% x0 = [alpha, beta, gamma, phih, mbarl, theta];
% load params;
% x0 = params(1:6);

A = [0,0,0,0,0,-1,1];
b = [1];
Aeq = [];
beq = [];
lb = [0.7,  1, 0.001, 1, 0.000001, 2,   2];
ub = [  1, 30,     1, 5,    0.5, 10, 5];

% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [0.7,  1, 0.001, 1, 0.000001, sigma];
% ub = [  1, 30,     1, 5,    0.5, 6];
nonlcon = @(x)mycon(x);
fun = @(x)calibration(x, data_moments);
% options = optimoptions('fmincon','Display','iter','MaxIterations', 50);
% [params,fval,exitflag,output] =  fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
options = optimoptions('patternsearch','Display','iter','MaxIterations', 50, 'MaxFunctionEvaluations', 300, 'Cache', 'on');
[params,fval,exitflag,output] =  patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
ptable = table(["alpha"; "beta"; "gamma"; "phih"; "mbar";"theta"; "sigma"], params');
ptable.Properties.VariableNames = ["Parameter","Value"];
ptable
save params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQUILIBRIUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load params;
[delta, lambda, s, g, Omega, iter, eq_error] = equilibrium(params);
model_moments = agmoments(delta, lambda, s, g, Omega, params);
mtable = table(["Ag Aquisition Prob"; "Av Scope"; "Relative Scope .1%"; "Rel Scale .1%"; "Ag Exit Rate"; "Growth Rate"; "Employment Growth"], [model_moments.aga, model_moments.avscope, model_moments.rel_scope01, model_moments.rel_scale01,  model_moments.agexit, model_moments.g, model_moments.empg]', data_moments');
mtable.Properties.VariableNames = ["Moment", "Model Value", "Data Value"];
mtable

% params(6) = 1.02;
% [delta2, lambda2, s2, g2, Omega2, iter, eq_error] = equilibrium(params);
% moments2 = agmoments(delta2, lambda2, s2, g2, Omega2, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUNTERFACTUALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params(8) = 0.98*params(8);
% [delta, lambda, s, g, Omega, iter, error] = equilibrium(guess, params);

function [c,ceq]  = mycon(x)
    c = x(1)*x(4)-0.999;
    ceq = [];
end