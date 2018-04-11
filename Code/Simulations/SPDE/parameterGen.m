function [p] = parameterGen(sigma,regime,noiseRegime,diffOn,singlePerc,stream)
    p = struct();
    
    baseExp = -3;
    
    p.prodLimit     = 200;
    
    if strcmp(noiseRegime,'MultIn')
        p.epsinMult= 2*62.5;
        p.epsinAdd   =0;
        p.epsoutMult =2*62.5;
        p.epsoutAdd  =0;
        p.epsRAdd    = 0;
        p.epsRMult =2*62.5;
        p.epsHAdd  =00;
        p.epsHMult =0;
        p.epsKAdd  =0;
        p.epsKMult =0;
    elseif strcmp(noiseRegime,'Add')
        p.epsinMult=0;
        p.epsinAdd   =6.25;
        p.epsoutMult =0;
        p.epsoutAdd  =2*6.25;
        p.epsRAdd    = 2*6.25;
        p.epsRMult =0;
        p.epsHAdd  =2*625;
        p.epsHMult =0;
        p.epsKAdd  =2*625;
        p.epsKMult =0;
    elseif strcmp(noiseRegime,'None')
        p.epsinMult=0;
        p.epsinAdd   =0;
        p.epsoutMult =0;
        p.epsoutAdd  =0;
        p.epsRAdd    = 0;
        p.epsRMult =0;
        p.epsHAdd  =00;
        p.epsHMult =0;
        p.epsKAdd  =0;
        p.epsKMult =0;
    elseif strcmp(noiseRegime,'AddSmall')
        p.epsinMult=0;
        p.epsinAdd   =20;
        p.epsoutMult =0;
        p.epsoutAdd  =0;
        p.epsRAdd    = 0;
        p.epsRMult =0;
        p.epsHAdd  =00;
        p.epsHMult =0;
        p.epsKAdd  =0;
        p.epsKMult =0;
    elseif strcmp(noiseRegime,'AddIn')
        p.epsinMult=0;
        p.epsinAdd   =.625;
        p.epsoutMult =0;
        p.epsoutAdd  =0;
        p.epsRAdd    = 0;
        p.epsRMult =0;
        p.epsHAdd  =00;
        p.epsHMult =0;
        p.epsKAdd  =0;
        p.epsKMult =0;
    end
    
    
    
	if diffOn
		p.DRA = 25.46;
	else
		p.DRA = 0;
	end
    if strcmp(regime,'Likun')
        p.Vbp    = (1e-6)/2*randOoM(sigma,stream);
        p.gamma  = 1e2*randOoM(sigma,stream);
        p.beta   = 1e1*randOoM(sigma,stream);
        p.kp     = 1e-4*randOoM(sigma,stream);
        p.kmax   = 1*randOoM(sigma,stream);
        p.beta0  = 400*randOoM(sigma,stream);
        p.lambda = 1e-1*randOoM(sigma,stream);
        p.kdeg   = (1e-2)*randOoM(sigma,stream);
        p.vmax   = (1e-5)*randOoM(sigma,stream);
        p.ron    = 3*randOoM(sigma,stream);
        p.roff   = 1e-3*randOoM(sigma,stream);
        p.mon    = 3*randOoM(sigma,stream);
        p.moff   = 1.3e-3*randOoM(sigma,stream);
        p.jalpha = 2*randOoM(sigma,stream);
        p.jbeta  = 1*randOoM(sigma,stream);
        p.bpdeg1 = 1e-4*randOoM(sigma,stream);
        p.bpdeg2 = 1e-4*randOoM(sigma,stream);
        p.Vr     = 1e-6*randOoM(sigma,stream);
        p.rdeg1  = 1e-4*randOoM(sigma,stream);
        p.rdeg2  = 1e-4*randOoM(sigma,stream);
		p.d      = 10*randOoM(sigma,stream);
		p.e      = 10*randOoM(sigma,stream);
    elseif strcmp(regime,'Flat')
        p.Vbp    = 1*randOoM(sigma,stream,baseExp);
        p.gamma  = 1*randOoM(sigma,stream,baseExp);
        p.beta   = 1*randOoM(sigma,stream,baseExp);
        p.kp     = 1*randOoM(sigma,stream,baseExp);
        p.kmax   = 1*randOoM(sigma,stream,baseExp);
        p.beta0  = 1*randOoM(sigma,stream,baseExp);
        p.lambda = 1*randOoM(sigma,stream,baseExp);
        p.kdeg   = 1*randOoM(sigma,stream,baseExp);
        p.vmax   = 1*randOoM(sigma,stream,baseExp); %boosted to stop negatives
        p.ron    = 1*randOoM(sigma,stream,baseExp);
        p.roff   = 1*randOoM(sigma,stream,baseExp);
        p.mon    = 1*randOoM(sigma,stream,baseExp);
        p.moff   = 1*randOoM(sigma,stream,baseExp);
        p.jalpha = 1*randOoM(sigma,stream,baseExp);
        p.jbeta  = 1*randOoM(sigma,stream,baseExp);
        p.bpdeg1 = 1*randOoM(sigma,stream,baseExp);
        p.bpdeg2 = 1*randOoM(sigma,stream,baseExp);
        p.Vr     = 1*randOoM(sigma,stream,baseExp);
        p.rdeg1  = 1*randOoM(sigma,stream,baseExp);
        p.rdeg2  = 1*randOoM(sigma,stream,baseExp);
        p.lambda = 1*randOoM(sigma,stream,baseExp);
		p.d      = 1*randOoM(sigma,stream,baseExp);
		p.e      = 1*randOoM(sigma,stream,baseExp);
		p.cH     = 1*randOoM(sigma,stream,baseExp);
		p.cK     = 1*randOoM(sigma,stream,baseExp);
		p.kappaH = 1*randOoM(sigma,stream,baseExp);
		p.kappaK = 1*randOoM(sigma,stream,baseExp);
		p.dH     = 1*randOoM(sigma,stream,baseExp);
		p.dK     = 1*randOoM(sigma,stream,baseExp);
		p.aH     = 1*randOoM(sigma,stream,baseExp);
		p.aK     = 1*randOoM(sigma,stream,baseExp);
        p.f0     = 500;
	elseif strcmp(regime,'FlatMod')
        p.Vbp    = 1*randOoM(sigma,stream,baseExp+3);
        p.gamma  = 1*randOoM(sigma,stream,baseExp);
        p.beta   = 1*randOoM(sigma,stream,baseExp+2);
        p.kp     = 1*randOoM(sigma,stream,baseExp);
        p.kmax   = 1*randOoM(sigma,stream,baseExp);
        p.beta0  = 1*randOoM(sigma,stream,baseExp);
        p.lambda = 1*randOoM(sigma,stream,baseExp);
        p.kdeg   = 1*randOoM(sigma,stream,baseExp+1);
        p.vmax   = 1*randOoM(sigma,stream,baseExp+1); %boosted to stop negatives
        p.ron    = 1*randOoM(sigma,stream,baseExp);
        p.roff   = 1*randOoM(sigma,stream,baseExp);
        p.mon    = 1*randOoM(sigma,stream,baseExp);
        p.moff   = 1*randOoM(sigma,stream,baseExp);
        p.jalpha = 1*randOoM(sigma,stream,baseExp);
        p.jbeta  = 1*randOoM(sigma,stream,baseExp);
        p.bpdeg1 = 1*randOoM(sigma,stream,baseExp+2);
        p.bpdeg2 = 1*randOoM(sigma,stream,baseExp);
        p.Vr     = 1*randOoM(sigma,stream,baseExp+3);
        p.rdeg1  = 1*randOoM(sigma,stream,baseExp+2);
        p.rdeg2  = 1*randOoM(sigma,stream,baseExp+2);
        p.lambda = 1*randOoM(sigma,stream,baseExp);
		p.d      = 1*randOoM(sigma,stream,baseExp);
		p.e      = 1*randOoM(sigma,stream,baseExp);
		p.cH     = 1*randOoM(sigma,stream,baseExp);
		p.cK     = 1*randOoM(sigma,stream,baseExp);
		p.kappaH = 1*randOoM(sigma,stream,baseExp);
		p.kappaK = 1*randOoM(sigma,stream,baseExp);
		p.dH     = 1*randOoM(sigma,stream,baseExp);
		p.dK     = 1*randOoM(sigma,stream,baseExp);
		p.aH     = 1*randOoM(sigma,stream,baseExp);
		p.aK     = 1*randOoM(sigma,stream,baseExp);
        p.f0     = 500;
    elseif strcmp(regime,'Mod')
        p.Vbp    = (1e-6)/2*randOoM(sigma,stream);
        p.gamma  = 1e2*randOoM(sigma,stream);
        p.beta   = 1e1*randOoM(sigma,stream);
        p.kp     = 1e-4*randOoM(sigma,stream);
        p.kmax   = 1*randOoM(sigma,stream);
        p.beta0  = 400*randOoM(sigma,stream);
        p.lambda = 1e-1*randOoM(sigma,stream);
        p.kdeg   = (1e-2)*randOoM(sigma,stream);
        p.vmax   = (1e-3)*randOoM(sigma,stream);
        p.ron    = 3*randOoM(sigma,stream);
        p.roff   = 1e-3*randOoM(sigma,stream);
        p.mon    = 3*randOoM(sigma,stream);
        p.moff   = 1.3e-3*randOoM(sigma,stream);
        p.jalpha = 2*randOoM(sigma,stream);
        p.jbeta  = 1*randOoM(sigma,stream);
        p.bpdeg1 = 1e-4*randOoM(sigma,stream);
        p.bpdeg2 = 1e-4*randOoM(sigma,stream);
        p.Vr     = 1e-6*randOoM(sigma,stream);
        p.rdeg1  = 1e-4*randOoM(sigma,stream);
        p.rdeg2  = 1e-4*randOoM(sigma,stream);
		p.d      = (1e-6)/2*randOoM(sigma,stream);
		p.e      = 1*randOoM(sigma,stream);
        p.cH     = 10*randOoM(sigma,stream);
		p.cK     = 10*randOoM(sigma,stream);
		p.kappaH = 10*randOoM(sigma,stream);
		p.kappaK = 10*randOoM(sigma,stream);
		p.dH     = 10*randOoM(sigma,stream);
		p.dK     = 10*randOoM(sigma,stream);
		p.aH     = 10*randOoM(sigma,stream);
		p.aK     = 10*randOoM(sigma,stream);
    end
    if singlePerc
        Names = fieldnames(p);
        %# Loop over the field-names and apply the function to each field
        for n = 1:length(Names)
            p.(Names{n}) =single(p.(Names{n}));
        end
    end
end