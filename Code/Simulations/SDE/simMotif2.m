function [RA,RAR,RAinit,RARinit,tRAv,tRARv] = simMotif2(vars,crabp,cypmax,N,dt,stream)
beta = vars(1);
delta= vars(2);
sigma= vars(3);
eta=   vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);
alpha = alpha0 * cypmax;
omega = omega0;
gamma = gamma0* crabp;
delta = delta * crabp;

RA = zeros(N,1);
RAR= zeros(N,1);

RAinit = (sqrt(2*beta*gamma*delta*omega*(2*alpha+eta)+beta*beta*gamma*gamma+ delta*delta*eta*eta*omega*omega) + beta*gamma + delta*eta*omega)/(2*gamma*(alpha+eta));
RARinit= (sqrt(4*beta*gamma*delta*omega*(alpha+eta) + (beta*gamma - delta * eta * omega)^2) + beta*gamma - delta*eta*omega) / (2*delta*(alpha+eta));

A = RAinit;
B = RARinit;
tRAv = (1/2).*(A.*gamma+delta.*omega).*(A.*gamma.*(alpha+delta+eta+gamma) ...
  +delta.*(delta+eta+gamma).*omega).^(-1).*(A.^2.*(alpha+eta).* ...
  gamma.^2+2.*A.*delta.*(alpha+eta).*gamma.*omega+delta.^2.*eta.* ...
  omega.^2).^(-1).*(A.^2.*(alpha+delta+eta).*gamma.^2+2.*A.*delta.*( ...
  alpha+delta+eta).*gamma.*omega+delta.^2.*(delta+eta).*omega.^2).* ...
  sigma.^2;

tRARv = (1/2).*gamma.^2.*(B+omega).^3.*(B.*(alpha+delta+eta+gamma)+(delta+ ...
  eta+gamma).*omega).^(-1).*(B.^2.*delta.*(alpha+eta)+B.*delta.*( ...
  alpha+2.*eta).*omega+omega.*(A.*alpha.*gamma+delta.*eta.*omega)) ...
  .^(-1).*sigma.^2;

RA(1) = RAinit;
RAR(1)= RARinit;

for i=2:N
    RA(i) = RA(i-1)  + (beta - (alpha*RAR(i-1)/(omega+RAR(i-1)) + gamma + eta )*RA(i-1) + delta*RAR(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RAR(i)= RAR(i-1) +(gamma*RA(i-1) - delta*RAR(i-1))*dt;
end

%{
fprintf('\nMean of RA is     %.2d \n',mean(RA))
fprintf('Variance of RA is %.2d',var(RA))
fprintf('\nVar/Mean of RA is %.2d \n',var(RA)/mean(RA))
%}
end