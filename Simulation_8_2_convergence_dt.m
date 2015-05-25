name = '8_2_conv_dt';
n = 32;
alpha = 10;
Rnucleus = 50;
V0 = -1;
x0 = -128;
sigmanorm = 0.05;
dt = 0.1;
tfinal = 30;
ndx = 1024;
nbsim = 100;
dt_max = 1;

Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 3, nbsim, 0, dt_max);

AnalyseSchroedinger_convergence_dt(name);