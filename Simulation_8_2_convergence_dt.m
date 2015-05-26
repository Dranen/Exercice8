name = '8_2_conv_dt';
n = 32;
alpha = 30;
Rnucleus = 50;
V0 = -1;
x0 = -128;
sigmanorm = 0.05;
dt = 0.05;
tfinal = 500;
ndx = 1024;
nbsim = 100;
dt_max = 0.5;

Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 1, 3, nbsim, 0, dt_max, 0);

AnalyseSchroedinger_convergence_dt(name);