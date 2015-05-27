%  fname='out_Delta100';
function AnalyseSchroedinger_convergence_dx(fname)
data=load([fname '_out.dat']);

dx=data(:,1);
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
loglog(dx,p,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('p')
saveas(gcf, [fname, '_dx_p.fig'])
saveas(gcf, [fname, '_dx_p.eps'], 'epsc')
figure
loglog(dx,p2,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('p^2')
saveas(gcf, [fname, '_dx_p2.fig'])
saveas(gcf, [fname, '_dx_p2.eps'], 'epsc')

figure
grid on
loglog(dx,errX,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('\Delta x')
saveas(gcf, [fname, '_dx_delta_x.fig'])
saveas(gcf, [fname, '_dx_delta_x.eps'], 'epsc')
figure
grid
loglog(dx,errP,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('\Delta p')
saveas(gcf, [fname, '_dx_delta_p.fig'])
saveas(gcf, [fname, '_dx_delta_p.eps'], 'epsc')

figure
grid on
loglog(dx,errX.*errP,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('(\Delta p)(\Delta x)')
saveas(gcf, [fname, '_dx_delta.fig'])
saveas(gcf, [fname, '_dx_delta.eps'], 'epsc')

figure
loglog(dx,Emean,'x')
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('E')
saveas(gcf, [fname, '_dx_E.fig'])
saveas(gcf, [fname, '_dx_E.eps'], 'epsc')

figure
loglog(dx,xmean,'k-','linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('<x>')
saveas(gcf, [fname, '_dx_x.fig'])
saveas(gcf, [fname, '_dx_x.eps'], 'epsc')

figure
loglog(dx,Pleft,'k-', dx,Pright,'r-', 'linewidth',lw);
set(gca,'fontsize',fs)
xlabel('dx')
ylabel('Probabilité')
saveas(gcf, [fname, '_dx_prob.fig'])
saveas(gcf, [fname, '_dx_prob.eps'], 'epsc')

%  xsup = -81.918;
%  fig=figure('Position', [0 00 641 300])
%  set(gca, 'fontsize', 12);
%  loglog(dx,abs((xmean-xsup)/xsup),'x');
%  xlabel('dx');
%  ylabel('erreur(<x>)');
%  fitobject = fit(log(dx),log(abs((xmean-xsup)/xsup)),'poly1');
%  hold on
%  loglog(dx,exp(feval(fitobject,log(dx))),'--r')
%  text(dx(round(end/2)),abs((xmean(round(end/2))-xsup)/xsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
%  set(gcf,'PaperPositionMode','auto')
%  
%  errXsup = -81.918;
%  fig=figure('Position', [0 00 641 300])
%  set(gca, 'fontsize', 12);
%  loglog(dx,abs((errX-errXsup)/errXsup),'x');
%  xlabel('dx');
%  ylabel('erreur(<x>)');
%  fitobject = fit(log(dx),log(abs((errX-errXsup)/errXsup)),'poly1');
%  hold on
%  loglog(dx,exp(feval(fitobject,log(dx))),'--r')
%  text(dx(round(end/2)),abs((errX(round(end/2))-errXsup)/errXsup),['pente = ',num2str(fitobject.p1)],'VerticalAlignment','Bottom')
%  set(gcf,'PaperPositionMode','auto')

 


