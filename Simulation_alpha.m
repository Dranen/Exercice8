name = '8_2_alpha';
n = 32;
Rnucleus = 50;
V0 = -1;
x0 = -128;
sigmanorm = 0.05;
dt = 0.5;
tfinal = 500;
ndx = 1024;
alpha = 0;

<<<<<<< HEAD
Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 4, 100, 0, 0, 1000);
=======
delete ('max_noy_alpha.dat');

for alpha=1:3:100
Simulation(name, n, alpha, Rnucleus, V0, x0, sigmanorm, dt, tfinal, ndx, 1, 1, 1, 0, 0);
alpha
end
>>>>>>> 767051b823762d95376b48ddf8556674237f7bb6

plot_alpha;