name = '8_2a';
n = 32;
Rnucleus = 50;
V0 = -1;
x0 = -128;
sigmanorm = 0.05;
dt = 0.5;
tfinal = 500;
ndx = 1024;

delete ('max_noy_alpha.dat');

for alpha=1:3:100
Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 1, 1, 0, 0);
alpha
end

plot_alpha;