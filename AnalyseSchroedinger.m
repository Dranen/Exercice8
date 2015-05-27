%  fname='out_Delta100';
function AnalyseSchroedinger(fname)
data=load([fname '_out.dat']);
fs=16; lw=1.5;
t=data(:,1);
Pleft=data(:,2);
Pright=data(:,3);
xmean=data(:,4);
xmean2=data(:,5);
p=data(:,6);
p2=data(:,7);
Emean=data(:,8);
errX=data(:,9);
errP=data(:,10);

figure
plot(t,p,'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t')
ylabel('p')
saveas(gcf, [fname, '_p.fig'])
saveas(gcf, [fname, '_p.eps'], 'epsc')
figure
plot(t(2:end),p2(2:end),'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t')
ylabel('p^2')
saveas(gcf, [fname, '_p2.fig'])
saveas(gcf, [fname, '_p2.eps'], 'epsc')

figure
grid on
plot(t,errX,'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t')
ylabel('\Delta x')
saveas(gcf, [fname, '_errX.fig'])
saveas(gcf, [fname, '_errX.eps'], 'epsc')
figure
grid
plot(t,errP,'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t')
ylabel('\Delta p')
saveas(gcf, [fname, '_errP.fig'])
saveas(gcf, [fname, '_errP.eps'], 'epsc')

figure
grid on
hold all
plot(t,errX.*errP,'linewidth',lw)
plot([min(t) max(t)], [0.5 0.5], 'k--')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('(\Delta p)(\Delta x)')
saveas(gcf, [fname, '_errPX.fig'])
saveas(gcf, [fname, '_errPX.eps'], 'epsc')

figure
plot(t,Emean,'linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t')
ylabel('E')
saveas(gcf, [fname, '_Emean.fig'])
saveas(gcf, [fname, '_Emean.eps'], 'epsc')

figure
m=1;
n=3;
p=0;
%  figure
subplot(m,n,p+1)
plot(t,Pleft,'k-', t,Pright,'r-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('t')
ylabel('P')

%  figure
subplot(m,n,p+2)
plot(t,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('x')
ylabel('<x>')

ErrorEnergyConserv=max(Emean)-min(Emean)

%-- 2D plot |psi(x,t)|^2
nt=length(t);
clear data
data=load([fname '_psi.dat']);
[nn,ii]=size(data);
nx=nn/nt;
xgrid=data(1:nx,2);
zpsiabs2=data(:,3);
psiabs2=reshape(zpsiabs2,nx,nt);
[X,T]=meshgrid(xgrid,t);
%  figure
subplot(m,n,p+3)
hs=contourf(X',T',psiabs2,20);
set(gca,'fontsize',fs)
xlabel('x')
ylabel('t')
shading flat

figure
contourf(X',T',psiabs2,20);
set(gca,'fontsize',fs)
colorbar
xlabel('x')
ylabel('t')
shading flat
saveas(gcf, [fname, '_psi.fig'])
saveas(gcf, [fname, '_psi.eps'], 'epsc')

figure
plot(t,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('x')
ylabel('<x>')
saveas(gcf, [fname, '_x.fig'])
saveas(gcf, [fname, '_x.eps'], 'epsc')

figure
plot(t,Pleft,'k-', t,Pright,'r-',t,Pleft+Pright,'b-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('t')
ylabel('Probabilité')
saveas(gcf, [fname, '_Proba.fig'])
saveas(gcf, [fname, '_Proba.eps'], 'epsc')

max(Pleft+Pright)
min(Pleft+Pright)

% figure
% for i = 1:max(size(t))
%     plot(X,psiabs2(:,i));
%     pause(0.01);
% end


