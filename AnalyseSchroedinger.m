%  fname='out_Delta100';
function AnalyseSchroedinger(fname)
data=load([fname '_out.dat']);

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
plot(t,p)
xlabel('t')
ylabel('p')
figure
plot(t(2:end),p2(2:end))
xlabel('t')
ylabel('p^2')

figure
grid on
plot(t,errX)
xlabel('t')
ylabel('\Delta x')
figure
grid
plot(t,errP)
xlabel('t')
ylabel('\Delta p')

figure
grid on
plot(t,errX.*errP)

figure
plot(t,Emean)
xlabel('t')
ylabel('E')

figure
fs=16; lw=1.5;
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

figure
plot(t,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('x')
ylabel('<x>')

figure
plot(t,Pleft,'k-', t,Pright,'r-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('t')
ylabel('P')

% figure
% for i = 1:max(size(t))
%     plot(X,psiabs2(:,i));
%     pause(0.01);
% end


