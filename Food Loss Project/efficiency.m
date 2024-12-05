rho = 0.03;
phi = 4;
gamma = 1;
pa = 100;
kappa = 10;
x = 0.5;
alpha = 0.5;
d0  = 1;
sigma = 0.05;
abar = 1;
z = 1.0100;
L  = 5;
beta = 1;
params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];

[l_eq, p_eq, iter, error] = equilibrium(params);
[moments] = agmoments(p_eq,l_eq,params);
w_eq = moments.W;
[w_sp, p_sp, l_sp] = planner(params);
[fl_fl, W_fl, p_fl, l_fl] = foodloss(params);

leq_array = zeros(1,11);   
peq_array = zeros(1,11);
Weq_array = zeros(1,11);
lsp_array = zeros(1,11);
psp_array = zeros(1,11);
Wsp_array = zeros(1,11);
lop_array = zeros(1,11);
pop_array = zeros(1,11);
Wop_array = zeros(1,11);
d_array = linspace(1,11, 11);
for i = 1:length(d_array)
    d0 = d_array(i);
    params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];
    [leq_array(i), peq_array(i), iter, error] = equilibrium(params);
    [moments] = agmoments(peq_array(i), leq_array(i),params);
    Weq_array(i) = moments.W;
    [Wsp_array(i), psp_array(i), lsp_array(i)] = planner(params);
    [flmin, Wop_array(i), pop_array(i), lop_array(i)] = foodloss(params);
end
leq_array = L - leq_array;
lsp_array = L - lsp_array;
lop_array = L - lop_array;
d0 = 1;
figure;
subplot(3,1,1);
plot(d_array, peq_array, d_array, psp_array, d_array, pop_array);
xlabel('Storage baseline');
ylabel('Price');
subplot(3,1,2);
plot(d_array, leq_array, d_array, lsp_array, d_array, lop_array);
xlabel('Storage baseline');
ylabel('Storage Investment');
subplot(3,1,3);
plot(d_array, Weq_array, d_array, Wsp_array, d_array, Wop_array);
xlabel('Storage baseline');
ylabel('Welfare');
saveas(gcf,'Output/SP_storage.png')

leq_array = zeros(1,11);
peq_array = zeros(1,11);
Weq_array = zeros(1,11);
lsp_array = zeros(1,11);
psp_array = zeros(1,11);
Wsp_array = zeros(1,11);
x_array = linspace(1,11, 11);
for i = 1:length(x_array)
    x = x_array(i);
    params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];
    [leq_array(i), peq_array(i), iter, error] = equilibrium(params);
    [moments] = agmoments(peq_array(i), leq_array(i),params);
    Weq_array(i) = moments.W;
    [Wsp_array(i), psp_array(i), lsp_array(i)] = planner(params);
end
leq_array = L - leq_array;
lsp_array = L - lsp_array;
x = 1;
figure;
subplot(3,1,1);
plot(x_array, peq_array, x_array, psp_array);
xlabel('Farm Size');
ylabel('Price');
subplot(3,1,2);
plot(x_array, leq_array, x_array, lsp_array);
xlabel('Farm Size');
ylabel('Storage Investment');
subplot(3,1,3);
plot(x_array, Weq_array, x_array, Wsp_array);
xlabel('Farm Size');
ylabel('Welfare');
saveas(gcf,'Output/SP_fsize.png')

leq_array = zeros(1,11);
peq_array = zeros(1,11);
Weq_array = zeros(1,11);
lsp_array = zeros(1,11);
psp_array = zeros(1,11);
Wsp_array = zeros(1,11);
z_array = linspace(3,23, 11);
for i = 1:length(z_array)
    z = z_array(i);
    params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];
    [leq_array(i), peq_array(i), iter, error] = equilibrium(params);
    [moments] = agmoments(peq_array(i), leq_array(i),params);
    Weq_array(i) = moments.W;
    [Wsp_array(i), psp_array(i), lsp_array(i)] = planner(params);
end
leq_array = L - leq_array;
lsp_array = L - lsp_array;
z = 10;
figure;
subplot(3,1,1);
plot(z_array, peq_array, z_array, psp_array);
xlabel('Home Production');
ylabel('Price');
subplot(3,1,2);
plot(z_array, leq_array, z_array, lsp_array);
xlabel('Home Production');
ylabel('Storage Investment');
subplot(3,1,3);
plot(z_array, Weq_array, z_array, Wsp_array);
xlabel('Home Production');
ylabel('Welfare');
saveas(gcf,'Output/SP_homeprod.png')

leq_array = zeros(1,11);
peq_array = zeros(1,11);
Weq_array = zeros(1,11);
lsp_array = zeros(1,11);
psp_array = zeros(1,11);
Wsp_array = zeros(1,11);
L_array = linspace(5,20, 11);
for i = 1:length(L_array)
    L = L_array(i);
    params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];
    [leq_array(i), peq_array(i), iter, error] = equilibrium(params);
    [moments] = agmoments(peq_array(i), leq_array(i),params);
    Weq_array(i) = moments.W;
    [Wsp_array(i), psp_array(i), lsp_array(i)] = planner(params);
end
leq_array = L_array - leq_array;
lsp_array = L_array - lsp_array;
L = 5;
figure;
subplot(3,1,1);
plot(L_array, peq_array, L_array, psp_array);
xlabel('Endownment');
ylabel('Price');
subplot(3,1,2);
plot(L_array, leq_array, L_array, lsp_array);
xlabel('Endownment');
ylabel('Storage Investment');
subplot(3,1,3);
plot(L_array, Weq_array, L_array, Wsp_array);
xlabel('Endownment');
ylabel('Welfare');
saveas(gcf,'Output/SP_homeprod.png')


leq_array = zeros(1,11);
peq_array = zeros(1,11);
Weq_array = zeros(1,11);
lsp_array = zeros(1,11);
psp_array = zeros(1,11);
Wsp_array = zeros(1,11);
a_array = linspace(0,1, 11);
for i = 1:length(L_array)
    abar = a_array(i);
    params = [rho, phi, gamma, pa, kappa, x, alpha, d0, sigma, z, L, abar, beta];
    [leq_array(i), peq_array(i), iter, error] = equilibrium(params);
    [moments] = agmoments(peq_array(i), leq_array(i),params);
    Weq_array(i) = moments.W;
    [Wsp_array(i), psp_array(i), lsp_array(i)] = planner(params);
end
leq_array = L - leq_array;
lsp_array = L - lsp_array;
abar = 1;
figure;
subplot(3,1,1);
plot(a_array, peq_array, a_array, psp_array);
xlabel('Mininum Consumption');
ylabel('Price');
subplot(3,1,2);
plot(a_array, leq_array, a_array, lsp_array);
xlabel('Mininum Consumption');
ylabel('Storage Investment');
subplot(3,1,3);
plot(a_array, Weq_array, a_array, Wsp_array);
xlabel('Mininum Consumption');
ylabel('Welfare');
saveas(gcf,'Output/SP_mincon.png')