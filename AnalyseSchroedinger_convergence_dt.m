%  fname='out_Delta100';
function AnalyseSchroedinger_convergence_dt(fname)
data=load([fname '_out.dat']);

dt=data(:,1);
Pleft=data(:,2);
Pright=data(:,3);
xmean=data(:,4);
xmean2=data(:,5);
p=data(:,6);
p2=data(:,7);
Emean=data(:,8);
errX=data(:,9);
errP=data(:,10);

fs=16; lw=1.5;
m=1;
n=3;

figure
loglog(dt,p)
xlabel('dt')
ylabel('p')
figure
loglog(dt,p2)
xlabel('dt')
ylabel('p^2')

figure
grid on
loglog(dt,errX)
xlabel('dt')
ylabel('\Delta x')
figure
grid
loglog(dt,errP)
xlabel('dt')
ylabel('\Delta p')

figure
grid on
loglog(dt,errX.*errP)

figure
loglog(dt,Emean)
xlabel('dt')
ylabel('E')

figure
loglog(dt,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('<x>')

figure
loglog(dt,Pleft,'k-', dt,Pright,'r-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('P')

 xsup = -81.918;
 fig=figure('Position', [0 00 641 300])
 set(gca, 'fontsize', 12);
 loglog(dt,abs((xmean-xsup)/xsup),'x');
 xlabel('dt');
 ylabel('erreur(<x>)');
 fitobject = fit(log(dt),log(abs((xmean-xsup)/xsup)),'poly1');
 hold on
 loglog(dt,exp(feval(fitobject,log(dt))),'--r')
 text(dt(round(end/2)),abs((xmean(round(end/2))-xsup)/xsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
 set(gcf,'PaperPositionMode','auto')
 
 errXsup = -81.918;
 fig=figure('Position', [0 00 641 300])
 set(gca, 'fontsize', 12);
 loglog(dt,abs((errX-errXsup)/errXsup),'x');
 xlabel('dt');
 ylabel('erreur(<x>)');
 fitobject = fit(log(dt),log(abs((errX-errXsup)/errXsup)),'poly1');
 hold on
 loglog(dt,exp(feval(fitobject,log(dt))),'--r')
 text(dt(round(end/2)),abs((errX(round(end/2))-errXsup)/errXsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
 set(gcf,'PaperPositionMode','auto')

 


