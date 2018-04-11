function [p] = convertVarsToP(vars)
	p = struct();

    p.vmax   = vars(1); %beta
    p.moff   = vars(2); %delta
    
    p.rdeg1  = vars(5); %eta
    p.kdeg   = vars(6); %alpha
    p.gamma  = vars(7); %omega
    p.mon    = vars(8); %gamma
    p.jalpha = vars(9); %nu
    p.jbeta  = vars(10); %lambda
    p.rdeg2  = vars(11); %r
    p.Vr     = vars(12); %zeta
    p.Vbp    = vars(13); %a
    p.bpdeg1 = vars(14); %u
    p.beta   = vars(15); %b
    p.kp     = vars(16); %c
    p.d        = vars(17); %d
    p.e        = vars(18); %e
    
    
    p.epsinMult = vars(3); %sigma
    p.epsRMult =vars(19); %sigma2
    
    
    p.bpdeg2 = 0;
    p.DRA = 200;
    
    p.epsinAdd = 0;
    p.epsoutMult = 0;
    p.epsoutAdd = 0;
    p.epsRAdd = 0;
    p.epsHAdd = 0;
    p.epsHMult = 0;
    p.epsKAdd = 0;
    p.epsKMult = 0;
    p.prodLimit = 200;
    p.cH = 2.7531e-07;
    p.cK = 3.0530e-06;
    p.kappaH = 9.9870e-07;
    p.kappaK = 1.0394e-07;
    p.dH = 3.7919e-06;
    p.dK = 1.4014e-07;
    p.aH = 3.2693e-07;
    p.aK = 9.1342e-06;
    p.f0 = 500;
    p.kmax=0;
    p.beta0=0;
    p.lambda=0;
    p.ron=0;
    p.roff=0;
end