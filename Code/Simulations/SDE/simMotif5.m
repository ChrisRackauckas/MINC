function [RA,RAR,RAinit,RARinit,tRAv,tRARv] = simMotif5(vars,crabp,cypmax,N,dt,stream)
% RA + BP upregulation
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
    d      =vars(17);
    e      =vars(18);
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
    RABPinit = (1/2).*delta.^(-1).*(alpha+eta).^(-1).*nu.^(-1).*u.^(-1).*zeta.^( ...
      -1).*((-1).*omega.*(delta.*e.*eta.*lambda.*r+beta.*gamma.*nu.* ...
      zeta)+e.*(delta.*e.*(alpha+eta).*lambda.*r+beta.*gamma.*nu.*zeta)) ...
      .^(-1).*((-1).*a.*((-1).*delta.*e.^2.*(alpha+eta).*lambda.*r+ ...
      beta.*gamma.*nu.*omega.*zeta+e.*(delta.*eta.*lambda.*omega.*r+(-1) ...
      .*beta.*gamma.*nu.*zeta)).*((-1).*delta.*eta.*lambda.*omega.*r+ ...
      beta.*gamma.*nu.*zeta+(delta.^2.*eta.^2.*lambda.^2.*omega.^2.* ...
      r.^2+2.*beta.*delta.*(2.*alpha+eta).*gamma.*lambda.*nu.*omega.*r.* ...
      zeta+beta.^2.*gamma.^2.*nu.^2.*zeta.^2).^(1/2))+d.*((-1).*beta.* ...
      gamma.*nu.*omega.*zeta.*((-1).*delta.*eta.*lambda.*omega.*r+beta.* ...
      gamma.*nu.*zeta+(delta.^2.*eta.^2.*lambda.^2.*omega.^2.*r.^2+2.* ...
      beta.*delta.*(2.*alpha+eta).*gamma.*lambda.*nu.*omega.*r.*zeta+ ...
      beta.^2.*gamma.^2.*nu.^2.*zeta.^2).^(1/2))+e.*(beta.^2.*gamma.^2.* ...
      nu.^2.*zeta.^2+delta.*eta.*lambda.*omega.*r.*(delta.*eta.*lambda.* ...
      omega.*r+(-1).*(delta.^2.*eta.^2.*lambda.^2.*omega.^2.*r.^2+2.* ...
      beta.*delta.*(2.*alpha+eta).*gamma.*lambda.*nu.*omega.*r.*zeta+ ...
      beta.^2.*gamma.^2.*nu.^2.*zeta.^2).^(1/2))+beta.*gamma.*nu.*zeta.* ...
      (2.*alpha.*delta.*lambda.*omega.*r+(delta.^2.*eta.^2.*lambda.^2.* ...
      omega.^2.*r.^2+2.*beta.*delta.*(2.*alpha+eta).*gamma.*lambda.*nu.* ...
      omega.*r.*zeta+beta.^2.*gamma.^2.*nu.^2.*zeta.^2).^(1/2)))));
    RARinit  = (gamma*zeta*nu)/(r*delta*lambda) * RAinit;
    Rinit    = zeta/r;
    BPinit   = ((-2).*omega.*u.*(delta.*e.*eta.*lambda.*r+beta.*gamma.*nu.*zeta)+ ...
      2.*e.*u.*(delta.*e.*(alpha+eta).*lambda.*r+beta.*gamma.*nu.*zeta)) ...
      .^(-1).*(delta.*e.*lambda.*(2.*a.*e.*(alpha+eta)+(-1).*(2.*a+d).* ...
      eta.*omega).*r+beta.*gamma.*nu.*((2.*a+d).*e+(-2).*(a+d).*omega).* ...
      zeta+d.*e.*(4.*beta.*delta.*(alpha+eta).*gamma.*lambda.*nu.* ...
      omega.*r.*zeta+((-1).*delta.*eta.*lambda.*omega.*r+beta.*gamma.* ...
      nu.*zeta).^2).^(1/2));
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
    RAR(i)  = RAR(i-1)  + (nu*RABP(i-1)*R(i-1) - lambda*RAR(i-1)*BP(i-1))*dt;
    R(i)    = R(i-1)    + (zeta - nu*RABP(i-1)*R(i-1) + lambda*RAR(i-1)*BP(i-1) - r*R(i-1))*dt; 
    BP(i)   = BP(i-1)   + (a-lambda*BP(i-1)*RAR(i-1) - gamma*BP(i-1)*RA(i-1) + (delta + nu*R(i-1))*RABP(i-1) - u*BP(i-1) + ((d*RAR(i-1))/(e+RAR(i-1))))*dt;
    end
end