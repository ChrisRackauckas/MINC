function [RA,RAR,RAinit,RARinit,tRAv,tRARv,D] = simMotif7(vars,crabp,cypmax,N,dt,stream)
%beta = 0;%vars(1)*10^-7;
delta= vars(2)* crabp;
sigma= vars(3);
%eta=   0;%vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);
nu    =vars(9)* crabp;
lambda0=vars(10);
alpha  = alpha0  * cypmax;
omega  = omega0;
gamma  = gamma0;
lambda = lambda0;
zeta = vars(12);
f    = vars(13);

p1 = vars(14)*1e-8; %Necessary
p2 = vars(15)*1e-4; %Necessary
p3 = 0;%vars(11)*1e-7; 

RA   = zeros(N,1);
RABP = zeros(N,1);
RAR  = zeros(N,1);
D    = zeros(N,1);

RAinit = gamma.^(-2).*nu.^(-1).*(f+p2).^(-1).*(nu.*p1.*(f+p2)+(lambda+p1).* ...
  p2.*zeta).^(-1).*(nu.*p1.*(f+p2)+(lambda+p1).*(delta.*f+p2.*( ...
  delta+zeta))).*(alpha.*nu.*p1.*(f+p2)+(-1).*gamma.*omega.*(nu.* ...
  p1.*(f+p2)+(lambda+p1).*p2.*zeta)+alpha.*(lambda+p1).*(delta.*f+ ...
  p2.*(delta+zeta)));

RABPinit = gamma.^(-1).*nu.^(-1).*(lambda+p1).*(nu.*p1.*(f+p2)+(lambda+p1).* ...
  p2.*zeta).^(-1).*(alpha.*nu.*p1.*(f+p2)+(-1).*gamma.*omega.*(nu.* ...
  p1.*(f+p2)+(lambda+p1).*p2.*zeta)+alpha.*(lambda+p1).*(delta.*f+ ...
  p2.*(delta+zeta)));
RARinit = (gamma.*nu.*p1.*(f+p2)+gamma.*(lambda+p1).*p2.*zeta).^(-1).*( ...
  alpha.*nu.*p1.*(f+p2)+(-1).*gamma.*omega.*(nu.*p1.*(f+p2)+(lambda+ ...
  p1).*p2.*zeta)+alpha.*(lambda+p1).*(delta.*f+p2.*(delta+zeta)));
Dinit = gamma.^(-1).*nu.^(-1).*(lambda+p1).*(f+p2).^(-1).*zeta.*(nu.*p1.*( ...
  f+p2)+(lambda+p1).*p2.*zeta).^(-1).*(alpha.*nu.*p1.*(f+p2)+(-1).* ...
  gamma.*omega.*(nu.*p1.*(f+p2)+(lambda+p1).*p2.*zeta)+alpha.*( ...
  lambda+p1).*(delta.*f+p2.*(delta+zeta)));

tRAv = 0;

tRARv = 0;

RA(1) = RAinit;
RABP(1) = RABPinit;
RAR(1)= RARinit;
D(1)  = Dinit;

for i=2:N
    RA(i)   = RA(i-1)   + (alpha*RA(i-1)/(omega + RAR(i-1)) - (gamma )*RA(i-1) + delta*RABP(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RABP(i) = RABP(i-1) + (gamma*RA(i-1) + lambda*RAR(i-1) + f*D(i-1) - (delta + nu + zeta + p3)*RABP(i-1))*dt;
    RAR(i)  = RAR(i-1)  + (nu*RABP(i-1) - (lambda+p1)*RAR(i-1))*dt;
    D(i)    = D(i-1)    + (zeta*RABP(i-1) - (f+p2)*D(i-1))*dt;
end

%{
tend = (N-1)*dt;
figure
subplot(4,1,1)
plot(0:dt:tend,RA)
subplot(4,1,2)
plot(0:dt:tend,RABP)
subplot(4,1,3)
plot(0:dt:tend,RAR)
subplot(4,1,4)
plot(0:dt:tend,D)
%}

end