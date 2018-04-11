dt = .01;
a = 1;
b = 2;
N = 1e1;
U = zeros(N,1);
U(1) = 2;
sigma = .01;

for i=1:N
    U(i+1) = U(i) + dt*(b - a*U(i));
end
%figure
%plot(U)
U(1) = U(end);
N = 1e8;
dt = 1e-5;
for i=1:N
    U(i+1) = U(i) + dt*(b - a*U(i)) + sqrt(dt)*sigma*randn;
end
tEnd = N*dt
[m,v] = cumMeanVar(U);
v(end)

figure
subplot(2,2,1)
hold on
title('X_t vs t')
xlabel('t')
ylabel('X_t')
plot(U)
hold off
epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
hgexport(gcf,'XtvsT',epsfig,'Format','eps');

subplot(2,2,2)
hold on
title('X_t vs t, Zoomed In')
xlabel('t')
ylabel('X_t')
plot(U(1:5e5))
hold off
epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
hgexport(gcf,'XtvsTShort',epsfig,'Format','eps');

subplot(2,2,3)
hold on
title('Variance of X_t vs t')
xlabel('t')
ylabel('Var[X_t]')
plot(v)
hold off
epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
hgexport(gcf,'VXtvsT',epsfig,'Format','eps');

subplot(2,2,4)
hold on
title('Variance of X_t vs t, Zoomed in')
xlabel('t')
ylabel('Var[X_t]')
plot(v(1:5e5))
hold off
epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
hgexport(gcf,'VXtvsTShort',epsfig,'Format','eps');

figure
hold on
title('Mean of X_t vs t')
xlabel('t')
ylabel('E[X_t]')
plot(m)
hold off
epsfig = hgexport('factorystyle');
epsfig.Format = 'eps';
hgexport(gcf,'EXtvsT',epsfig,'Format','eps');
