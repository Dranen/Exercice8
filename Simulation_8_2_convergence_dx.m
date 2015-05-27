name = '8_2_conv_dx';
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
ndx_max = 4096;

Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 1, 3, nbsim, ndx_max, 0, 0);

AnalyseSchroedinger_convergence_dx(name);