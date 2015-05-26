name = '8_2a';
n = 32;
alpha = 30;
Rnucleus = 50;
V0 = -1;
x0 = -128;
sigmanorm = 0.05;
dt = 0.5;
tfinal = 500;
ndx = 1024;
echt = 2;

Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, echt, 1, 1, 1, 0, 0, 0);

AnalyseSchroedinger(name);