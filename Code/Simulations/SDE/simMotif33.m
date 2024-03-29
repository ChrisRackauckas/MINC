function [RA,RAR,RAinit,RARinit,tRAv,tRARv,RABP] = simMotif33(vars,crabp,cypmax,N,dt,stream)
beta = vars(1);
delta= vars(2);
sigma= vars(3);
eta=   vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);%*1e5; %high gamma high nu works, low gamma low nu works
nu    =vars(9);%*1e5;%*1e5;
lambda0=vars(10);
alpha  = alpha0  * cypmax;
omega  = omega0;
gamma  = gamma0  * crabp;
lambda = lambda0 * crabp;
sigma2 =vars(19);
sigma3 =vars(20);
RA   = zeros(N,1);
RABP = zeros(N,1);
RAR  = zeros(N,1);

RAinit   = (beta*gamma*nu - delta*eta*lambda*omega + sqrt(4*beta*gamma*delta*(alpha+eta)*lambda*nu*omega+(beta*gamma*nu-delta*eta*lambda*omega)^2))/(2*gamma*(alpha+eta)*nu);
RABPinit = (gamma/delta)*RAinit;
RARinit  = (nu/lambda)*RABPinit;

tRAv = (1/2).*delta.^3.*(lambda.*omega+(1/2).*delta.^(-1).*(alpha+eta).^( ...
  -1).*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.* ...
  delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1) ...
  .*delta.*eta.*lambda.*omega).^2).^(1/2))).*(delta.^2.*eta.* ...
  lambda.^2.*omega.^2+delta.*lambda.*omega.*(beta.*gamma.*nu+(-1).* ...
  delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).*gamma.* ...
  lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega).^2).^(1/2))+(1/4).*(alpha+eta).^(-1).*(beta.*gamma.*nu+(-1) ...
  .*delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).*gamma.* ...
  lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega).^2).^(1/2)).^2).^(-1).*(delta.^2.*lambda.^2.*(gamma.^2.*( ...
  lambda+nu)+(delta+lambda+nu).*((delta+eta).*(eta+lambda)+eta.*nu)+ ...
  gamma.*((lambda+nu).*(2.*eta+lambda+nu)+delta.*(eta+2.*lambda+nu)) ...
  ).*omega.^2+(1/2).*delta.*(alpha+eta).^(-1).*lambda.*(alpha.*( ...
  delta.^2+(lambda+nu).*(2.*eta+lambda+nu)+delta.*(2.*eta+lambda+2.* ...
  nu)+gamma.*(delta+2.*(lambda+nu)))+2.*(gamma.^2.*(lambda+nu)+( ...
  delta+lambda+nu).*((delta+eta).*(eta+lambda)+eta.*nu)+gamma.*(( ...
  lambda+nu).*(2.*eta+lambda+nu)+delta.*(eta+2.*lambda+nu)))).* ...
  omega.*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.* ...
  delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1) ...
  .*delta.*eta.*lambda.*omega).^2).^(1/2))+(1/4).*(alpha+eta).^(-2) ...
  .*(gamma.^2.*(lambda+nu)+alpha.^2.*(delta+lambda+nu)+(delta+ ...
  lambda+nu).*((delta+eta).*(eta+lambda)+eta.*nu)+gamma.*((lambda+ ...
  nu).*(2.*eta+lambda+nu)+delta.*(eta+2.*lambda+nu))+alpha.*((delta+ ...
  lambda+nu).*(delta+2.*eta+lambda+nu)+gamma.*(delta+2.*(lambda+nu)) ...
  )).*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.* ...
  delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1) ...
  .*delta.*eta.*lambda.*omega).^2).^(1/2)).^2).^(-1).*(delta.* ...
  lambda.^3.*(delta.^2.*(eta+lambda)+eta.*(lambda+nu).*(eta+gamma+ ...
  lambda+nu)+delta.*(gamma.*lambda+(eta+lambda).^2+(2.*eta+lambda).* ...
  nu)).*omega.^3+(1/2).*(alpha+eta).^(-1).*lambda.^2.*(3.*( ...
  delta.^2.*(eta+lambda)+eta.*(lambda+nu).*(eta+gamma+lambda+nu)+ ...
  delta.*(gamma.*lambda+(eta+lambda).^2+(2.*eta+lambda).*nu))+ ...
  alpha.*(2.*delta.^2+3.*delta.*(eta+lambda)+4.*delta.*nu+(lambda+ ...
  nu).*(3.*eta+2.*gamma+2.*(lambda+nu)))).*omega.^2.*(beta.*gamma.* ...
  nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).* ...
  gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.* ...
  lambda.*omega).^2).^(1/2))+(1/4).*delta.^(-1).*(alpha+eta).^(-2).* ...
  lambda.*(2.*alpha.^2.*(delta+lambda+nu)+3.*(delta.^2.*(eta+lambda) ...
  +eta.*(lambda+nu).*(eta+gamma+lambda+nu)+delta.*(gamma.*lambda+( ...
  eta+lambda).^2+(2.*eta+lambda).*nu))+alpha.*(3.*delta.^2+5.* ...
  delta.*(eta+lambda)+6.*delta.*nu+(lambda+nu).*(5.*eta+3.*gamma+3.* ...
  (lambda+nu)))).*omega.*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega+(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+( ...
  beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega).^2).^(1/2)).^2+( ...
  1/8).*delta.^(-2).*(alpha+eta).^(-3).*((alpha+delta+eta).*(gamma.* ...
  lambda+alpha.*(delta+lambda)+(delta+lambda).*(eta+lambda))+(( ...
  alpha+eta).*(alpha+2.*delta+eta+gamma)+(2.*alpha+delta+2.*eta).* ...
  lambda).*nu+(alpha+eta).*nu.^2).*(beta.*gamma.*nu+(-1).*delta.* ...
  eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.* ...
  nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega).^2).^( ...
  1/2)).^3).*sigma.^2;

tRARv = (1/2).*gamma.^2.*lambda.^(-1).*nu.^2.*(delta.*lambda.*omega+(1/2) ...
  .*(alpha+eta).^(-1).*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega+(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+( ...
  beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega).^2).^(1/2))).^3.* ...
  (lambda.*(delta+eta+gamma+lambda+nu).*omega+(1/2).*delta.^(-1).*( ...
  alpha+eta).^(-1).*(alpha+delta+eta+gamma+lambda+nu).*(beta.* ...
  gamma.*nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+ ...
  eta).*gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.* ...
  lambda.*omega).^2).^(1/2))).*(delta.^2.*eta.*lambda.^2.*omega.^2+ ...
  delta.*lambda.*omega.*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega+(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+( ...
  beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega).^2).^(1/2))+(1/4) ...
  .*(alpha+eta).^(-1).*(beta.*gamma.*nu+(-1).*delta.*eta.*lambda.* ...
  omega+(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.*nu.*omega+( ...
  beta.*gamma.*nu+(-1).*delta.*eta.*lambda.*omega).^2).^(1/2)).^2) ...
  .^(-1).*(delta.^2.*lambda.^2.*(gamma.^2.*(lambda+nu)+(delta+ ...
  lambda+nu).*((delta+eta).*(eta+lambda)+eta.*nu)+gamma.*((lambda+ ...
  nu).*(2.*eta+lambda+nu)+delta.*(eta+2.*lambda+nu))).*omega.^2+( ...
  1/2).*delta.*(alpha+eta).^(-1).*lambda.*(alpha.*(delta.^2+(lambda+ ...
  nu).*(2.*eta+lambda+nu)+delta.*(2.*eta+lambda+2.*nu)+gamma.*( ...
  delta+2.*(lambda+nu)))+2.*(gamma.^2.*(lambda+nu)+(delta+lambda+nu) ...
  .*((delta+eta).*(eta+lambda)+eta.*nu)+gamma.*((lambda+nu).*(2.* ...
  eta+lambda+nu)+delta.*(eta+2.*lambda+nu)))).*omega.*(beta.*gamma.* ...
  nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).* ...
  gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.* ...
  lambda.*omega).^2).^(1/2))+(1/4).*(alpha+eta).^(-2).*(gamma.^2.*( ...
  lambda+nu)+alpha.^2.*(delta+lambda+nu)+(delta+lambda+nu).*((delta+ ...
  eta).*(eta+lambda)+eta.*nu)+gamma.*((lambda+nu).*(2.*eta+lambda+ ...
  nu)+delta.*(eta+2.*lambda+nu))+alpha.*((delta+lambda+nu).*(delta+ ...
  2.*eta+lambda+nu)+gamma.*(delta+2.*(lambda+nu)))).*(beta.*gamma.* ...
  nu+(-1).*delta.*eta.*lambda.*omega+(4.*beta.*delta.*(alpha+eta).* ...
  gamma.*lambda.*nu.*omega+(beta.*gamma.*nu+(-1).*delta.*eta.* ...
  lambda.*omega).^2).^(1/2)).^2).^(-1).*sigma.^2;

RA(1) = RAinit;
RABP(1) = RABPinit;
RAR(1)= RARinit;

for i=2:N
    RA(i)   = RA(i-1)   + (beta - (alpha*RAR(i-1)/(omega+RAR(i-1)) + gamma + eta )*RA(i-1) + delta*RABP(i-1))*dt + sigma*sqrt(dt)*RA(i-1)*randn(stream);
    RABP(i) = RABP(i-1) + (gamma*RA(i-1) + lambda*RAR(i-1) - (delta + nu)*RABP(i-1))*dt;% + sigma3*sqrt(dt)*RABP(i-1)*randn(stream);
    RAR(i)  = RAR(i-1)  + (nu*RABP(i-1) - lambda*RAR(i-1))*dt + sigma2*sqrt(dt)*RAR(i-1)*randn(stream);
end

end

%{
tend = N*dt
figure
subplot(3,1,1)
plot(0:dt:tend-dt,RA)
subplot(3,1,2)
plot(0:dt:tend-dt,RABP)
subplot(3,1,3)
plot(0:dt:tend-dt,RAR)
%}

%{
fprintf('\nMean of RA is     %.2d \n',mean(RA))
fprintf('Variance of RA is %.2d \n',var(RA))
%fprintf('\nVar/Mean of RA is %.2d \n',var(RA)/mean(RA))
%}
