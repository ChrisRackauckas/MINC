function [RA,RAR,RAinit,RARinit,tRAv,tRARv] = simMotif1(vars,crabp,cypmax,N,dt,stream)
beta = vars(1);
delta= vars(2);
sigma= vars(3);
eta=   vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);
alpha = alpha0 * cypmax;
omega = omega0;
gamma = gamma0 * crabp;
delta = delta * crabp;

RA = zeros(N,1);
RAR= zeros(N,1);

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

%{
fprintf('\nMean of RA is     %.2d \n',mean(RA))
fprintf('Variance of RA is %.2d',var(RA))
fprintf('\nVar/Mean of RA is %.2d \n',var(RA)/mean(RA))
%}
end