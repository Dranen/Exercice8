

bon=load('max_noy_alpha.dat');

figure;
plot(bon(:,1),bon(:,2),'+');
hold on;
set(gca,'fontsize',20);
xlabel('alpha');
ylabel('max \psi^2 \in [R_{nucleus},0]');