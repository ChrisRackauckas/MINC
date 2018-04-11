dt = 1e-2;
N = 5e6;

stream = RandStream('mrg32k3a');

crabp = 1;

beta = 0.001;
delta= .1;
sigma= 3e-3;
eta=   0.0001;
gamma0=.1;
gamma = gamma0* crabp;
delta = delta * crabp;

RA = zeros(N,1);
RAR= zeros(N,1);
t  = zeros(N,1);

RAinit = beta / eta;
RARinit= beta * gamma / (delta * eta);

A = RAinit;
B = RARinit;
tRAv = (delta + eta)*sigma*sigma / (2*eta*(gamma + delta + eta));

tRARv = (gamma*gamma*sigma*sigma)/(2*delta*eta*(gamma+delta+eta));

RA(1) = RAinit;
RAR(1)= RARinit;

for i=2:N
    RA(i) = RA(i-1)  + (beta - (gamma + eta )*RA(i-1) + delta*RAR(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RAR(i)= RAR(i-1) +(gamma*RA(i-1) - delta*RAR(i-1))*dt;
end

crabp = .005;

gamma = gamma0* crabp;
delta = delta * crabp;

RA2 = zeros(N,1);
RAR2= zeros(N,1);

RA2(1) = RAinit;
RAR2(1)= RARinit;

for i=2:N
    RA2(i) = RA2(i-1)  + (beta - (gamma + eta )*RA2(i-1) + delta*RAR2(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RAR2(i)= RAR2(i-1) +(gamma*RA(i-1) - delta*RAR2(i-1))*dt;
    t(i) = t(i-1)+dt;
end

cov1 = cov(RA,RAR);
cov2 = cov(RA2,RAR2);

cov1 = cov1(1,2);
cov2 = cov2(1,2);


fntSze = 18;
lw1 = 1.5;
lw2 = 1.8;
figDefs = get(0,'defaultfigureposition');
margvHist = .1;
marghHist = .1;
margvHistBot = .065;

figure('Position',[figDefs(1),figDefs(2),720,480])
subplot_tight(2,1,1,[margvHist marghHist]);
hold on 
plot(t,RA,'linewidth',lw1)
plot(t,RAR,'linewidth',lw1)
title('(I) RA Timeseries','FontSize',fntSze)
legend('RA','RAR','location','southeast')
axis([0 N*dt 9.6 10.4])
set(gca,'XTickLabel',[])
set(gca,'FontSize',fntSze)
hold off
subplot_tight(2,1,2,[margvHist marghHist]);
hold on 
plot(t,RA2,'linewidth',lw1)
plot(t,RAR2,'linewidth',lw1)
title('(II) Decreased Binding Rates','FontSize',fntSze)
y = ylabel('Concentration (uM)','FontSize',fntSze);
x = xlabel('Time (s)','FontSize',fntSze);
legend('RA','RAR','location','southeast')
axis([0 N*dt 9.6 10.4])
set(gca,'FontSize',fntSze)
hold off

figure('Position',[figDefs(1),figDefs(2),400,780])
bar([1 2],[cov1 cov2])
set(gca,'FontSize',fntSze)
str = {'(I)','(II)'};
y = ylabel('Covariance','FontSize',fntSze);
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))