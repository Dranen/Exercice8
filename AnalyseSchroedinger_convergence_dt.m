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
loglog(dt,p,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('p')
saveas(gcf, [fname, '_dt_p.fig'])
saveas(gcf, [fname, '_dt_p.eps'], 'epsc')
figure
loglog(dt,p2,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('p^2')
saveas(gcf, [fname, '_dt_p2.fig'])
saveas(gcf, [fname, '_dt_p2.eps'], 'epsc')

figure
grid on
loglog(dt,errX,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('\Delta x')
saveas(gcf, [fname, '_dt_delta_x.fig'])
saveas(gcf, [fname, '_dt_delta_x.eps'], 'epsc')
figure
grid
loglog(dt,errP,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('\Delta p')
saveas(gcf, [fname, '_dt_delta_p.fig'])
saveas(gcf, [fname, '_dt_delta_p.eps'], 'epsc')

figure
grid on
loglog(dt,errX.*errP,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('(\Delta p)(\Delta x)')
saveas(gcf, [fname, '_dt_delta.fig'])
saveas(gcf, [fname, '_dt_delta.eps'], 'epsc')

figure
loglog(dt,Emean,'x')
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('E')
saveas(gcf, [fname, '_dt_E.fig'])
saveas(gcf, [fname, '_dt_E.eps'], 'epsc')

figure
loglog(dt,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('<x>')
saveas(gcf, [fname, '_dt_x.fig'])
saveas(gcf, [fname, '_dt_x.eps'], 'epsc')

figure
loglog(dt,Pleft,'k-', dt,Pright,'r-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dt')
ylabel('Probabilité')
saveas(gcf, [fname, '_dt_prob.fig'])
saveas(gcf, [fname, '_dt_prob.eps'], 'epsc')

%  xsup = -81.918;
%  fig=figure('Position', [0 00 641 300])
%  set(gca, 'fontsize', 12);
%  loglog(dt,abs((xmean-xsup)/xsup),'x');
%  xlabel('dt');
%  ylabel('erreur(<x>)');
%  fitobject = fit(log(dt),log(abs((xmean-xsup)/xsup)),'poly1');
%  hold on
%  loglog(dt,exp(feval(fitobject,log(dt))),'--r')
%  text(dt(round(end/2)),abs((xmean(round(end/2))-xsup)/xsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
%  set(gcf,'PaperPositionMode','auto')
%  
%  errXsup = -81.918;
%  fig=figure('Position', [0 00 641 300])
%  set(gca, 'fontsize', 12);
%  loglog(dt,abs((errX-errXsup)/errXsup),'x');
%  xlabel('dt');
%  ylabel('erreur(<x>)');
%  fitobject = fit(log(dt),log(abs((errX-errXsup)/errXsup)),'poly1');
%  hold on
%  loglog(dt,exp(feval(fitobject,log(dt))),'--r')
%  text(dt(round(end/2)),abs((errX(round(end/2))-errXsup)/errXsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
%  set(gcf,'PaperPositionMode','auto')

 


