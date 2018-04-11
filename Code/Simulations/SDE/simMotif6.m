function [RA,RAR,RAinit,RARinit,tRAv,tRARv] = simMotif6(vars,crabp,cypmax,N,dt,stream)
beta = vars(1);
delta= vars(2);
sigma= vars(3);
eta=   vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);
nu    =vars(9);
lambda0=vars(10);
alpha  = alpha0  * cypmax;
omega  = omega0;
gamma  = gamma0  * crabp;
lambda = lambda0 * crabp;

RA   = zeros(N,1);
RABP = zeros(N,1);
RAR  = zeros(N,1);

RAinit   = (1/2).*eta.^(-1).*gamma.^(-1).*nu.^(-1).*(alpha.*delta.*lambda+ ...
  beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.*delta.* ...
  eta.*gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+delta.*lambda.*( ...
  alpha+(-1).*eta.*omega)).^2).^(1/2));
RABPinit = (gamma/delta)*RAinit;
RARinit  = (nu/lambda)*RABPinit;

tRAv = 0;

tRARv = 0;

RA(1) = RAinit;
RABP(1) = RABPinit;
RAR(1)= RARinit;

for i=2:N
    RA(i)   = RA(i-1)   + (beta + ((alpha*RA(i-1))/(omega+RAR(i-1))) - ( gamma + eta )*RA(i-1) + delta*RABP(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RABP(i) = RABP(i-1) + (gamma*RA(i-1) + lambda*RAR(i-1) - (delta + nu)*RABP(i-1))*dt;
    RAR(i)  = RAR(i-1)  + (nu*RABP(i-1) - lambda*RAR(i-1))*dt;
end

end

%{
tend = (N-1)*dt;
figure
subplot(3,1,1)
plot(0:dt:tend,RA)
subplot(3,1,2)
plot(0:dt:tend,RABP)
subplot(3,1,3)
plot(0:dt:tend,RAR)
%}

%{
fprintf('\nMean of RA is     %.2d \n',mean(RA))
fprintf('Variance of RA is %.2d \n',var(RA))
%fprintf('\nVar/Mean of RA is %.2d \n',var(RA)/mean(RA))
%}