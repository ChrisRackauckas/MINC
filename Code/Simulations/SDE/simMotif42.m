function [RA,RAR,RAinit,RARinit,tRAv,tRARv] = simMotif42(vars,crabp,cypmax,N,dt,stream)
%Motif 4 with additive noise on RA, multiplicative noise on RAR
beta = vars(1);
delta= vars(2);
sigma= vars(3);
eta=   vars(5);
alpha0=vars(6);
omega0=vars(7);
gamma0=vars(8);
nu    =vars(9);
lambda0=vars(10);
r      =vars(11);
zeta   =vars(12);
a      =vars(13);
u      =vars(14);
b      =vars(15);
c      =vars(16);
sigma2 =vars(19);
alpha  = alpha0  * cypmax;
omega  = omega0;
gamma  = gamma0;
lambda = lambda0;
a      = a*crabp;

RA   = zeros(N,1);
RABP = zeros(N,1);
RAR  = zeros(N,1);
R    = zeros(N,1);
BP   = zeros(N,1);
RAout= zeros(N,1);

RAinit   = (beta*gamma*zeta*nu - r*delta*eta*lambda*omega + sqrt(4*r*beta*gamma*delta*zeta*(alpha+eta)*lambda*nu*omega + (beta*gamma*zeta*nu - r*delta*eta*lambda*omega)^2))/...
    (2*gamma*zeta*(alpha+eta)*nu);
RABPinit = ((a*gamma)/(u*delta)) * RAinit;
RARinit  = ((u*zeta*nu)/(a*r*lambda)) * RABPinit;
Rinit    = zeta/r;
BPinit   = a/u;
RAoutinit= (beta + c*RAinit)/b;

RA(1)   = RAinit;
RABP(1) = RABPinit;
RAR(1)  = RARinit;
R(1)    = Rinit;
BP(1)   = BPinit;
RAout(1)=RAoutinit;

tRAv = 0;
tRARv = 0;

for i=2:N
	RAout(i)= RAout(i-1)+ (beta - b*RAout(i-1) + c*RA(i-1))*dt;
    RA(i)   = RA(i-1)   + (b*RAout(i-1) - c*RA(i-1) - (alpha*RAR(i-1)/(omega+RAR(i-1)) + gamma*BP(i-1) + eta )*RA(i-1) + delta*RABP(i-1))*dt + sigma*sqrt(dt)*randn(stream);
    RABP(i) = RABP(i-1) + (gamma*RA(i-1)*BP(i-1) + lambda*RAR(i-1)*BP(i-1) - (delta + nu*R(i-1))*RABP(i-1))*dt;
    RAR(i)  = RAR(i-1)  + (nu*RABP(i-1)*R(i-1) - lambda*RAR(i-1)*BP(i-1))*dt + sigma2*sqrt(dt)*RAR(i-1)*randn(stream);
    R(i)    = R(i-1)    + (zeta - nu*RABP(i-1)*R(i-1) + lambda*RAR(i-1)*BP(i-1) - r*R(i-1))*dt;
    BP(i)   = BP(i-1)   + (a-lambda*BP(i-1)*RAR(i-1) - gamma*BP(i-1)*RA(i-1) + (delta + nu*R(i-1))*RABP(i-1) - u*BP(i-1))*dt;
end

end
