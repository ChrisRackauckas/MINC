function [RAinsave,RARsave,sharpData,gradSharpData,RAoutsave,Rsave,BPsave,RABPsave,RAinGrad,RAoutGrad,RARGrad,RGrad,BPGrad,RABPGrad,HoxGrad,KroxGrad] = main(tEnd,dx,dy,dt,bpMult,cypMult,p,reaction,initRun,tTot2,dt2,stream,gpuID,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,gradSharpData,sharpLocations,useInits,knock)

    %%% Parameters %%%
    xStart = -100;
    xEnd   = 400;
    yEnd   = 50;
    
    %Calculuate Some Numbers
    m = (xEnd-xStart) / dx + 1; %Starts at 0
    n = yEnd / dy +1;
    tTot = tEnd / dt + 1;
    X = repmat(xStart:dx:xEnd,n,1);
    Vbp    = bpMult * cypMult^vbpAdjust * p.Vbp;
    gamma  = p.gamma;
    beta   = p.beta;
    kp     = p.kp;
    kmax   = p.kmax;
    beta0  = p.beta0;
    lambda = p.lambda;
    prodLimit = p.prodLimit;
    kdeg   = cypMult * p.kdeg;
    vmax   = p.vmax;
    ron    = p.ron;
    roff   = p.roff;
    mon    = p.mon;
    moff   = p.moff;
    jalpha = p.jalpha;
    jbeta  = p.jbeta;
    bpdeg1 = p.bpdeg1;
    bpdeg2 = p.bpdeg2;
    Vr     = p.Vr;
    rdeg1  = p.rdeg1;
    rdeg2  = p.rdeg2;   
	d      = p.d;
	e      = p.e;
    cH     = p.cH;
    cK     = p.cK;
    kappaH = p.kappaH;
    kappaK = p.kappaK;
    dH     = p.dH;
    dK     = p.dK;
    DRA    = p.DRA;
    lambda = p.lambda;
    f0     = p.f0;
	
	epsinMult = p.epsinMult ;
    epsinAdd  = p.epsinAdd  ;
    epsoutMult= p.epsoutMult;
    epsoutAdd = p.epsoutAdd ;
	epsRAdd   = p.epsRAdd   ;
    epsRMult  = p.epsRMult  ;
    epsHAdd   = p.epsHAdd   ;
    epsHMult  = p.epsHMult  ;
    epsKAdd   = p.epsKAdd   ;
    epsKMult  = p.epsKMult  ;
	
    

    %%% Script %%%
    production = vmax*(X >= prodLimit);
    alphaX = DRA/(dx*dx);
    alphaY = DRA/(dy*dy);
    alpha  = sqrt(dt) / (dx*dx*dy*dy);
    checkX1 = round(1*m/5);
    checkX2 = round(2*m/5);
    checkY  = round(n/2);
	checkX3 = round(3*m/5);
	checkX4 = round(4*m/5);
	gpuEnabled = gpuID>=1;
	sharpData = zeros(length(sharpLocations)*6,1);
    
    %Setup Initial Value Matrix
    highlim = 255;
    lowlim  = 105;
    interp = (X-lowlim)/(highlim-lowlim).*(X <highlim & X>=lowlim);
    
    if useInits
        inits = getInits(p,knock);
    end
    
    if bpMult ==1
        RAinit = inits.RAin;
        RAoutinit = inits.RAout;
        Rinit = inits.R;
        BPinit = inits.BP;
        RARinit = inits.RAR;
        RABPinit = inits.RABP;
        Rlowinit = inits.Rlow;
        BPlowinit = inits.BPlow;
    else
        RAinit = inits.RAin2;
        RAoutinit = inits.RAout2;
        Rinit = inits.R2;
        BPinit = inits.BP2;
        RARinit = inits.RAR2;
        RABPinit = inits.RABP2;
        Rlowinit = inits.Rlow2;
        BPlowinit = inits.BPlow2;
    end
    
    if DRA == 0
        RAout = RAoutinit*(X >= prodLimit);
        RAin  = RAinit*(X >= prodLimit) ;
        R     = Rlowinit*(X < prodLimit) + Rinit*(X >= prodLimit);
        RAR   = RARinit*(X >= prodLimit) ;
        BP    = BPlowinit*(X < prodLimit) + BPinit*(X >= prodLimit);
        RABP  = RABPinit*(X >= prodLimit) ;
    else
        RAout = RAoutinit*(X >= highlim) + RAoutinit.*interp;
        RAin  = RAinit*(X >= highlim) + RAinit.*interp;
        R     = Rlowinit*(X < highlim) + Rinit*(X >= prodLimit) + (Rinit.*interp + Rlowinit.*(1-interp)).*(X <highlim & X>=lowlim);
        RAR   = RARinit*(X >= highlim) + RARinit.*interp;
        BP    = BPlowinit*(X < highlim) + BPinit*(X >= prodLimit) + (BPinit.*interp + BPlowinit.*(1-interp)).*(X <highlim & X>=lowlim);
        RABP  = RABPinit*(X >= highlim) + RABPinit.*interp;
    end
        
        
    Hox   = .1*rand(n,m);
    Krox  = zeros(n,m);

    RAoutNoise = 0;
    RAinNoise  = 0;
    RARNoise   = 0;
    
    
    printSpacing = 10000;
    
    %Detailed Saves
    RAinsave = zeros(tTot2,4);
    RARsave  = zeros(tTot2,4);
    
    
    if initRun 
        RAoutsave= zeros(tTot2,1);
        Rsave    = zeros(tTot2,1);
        BPsave   = zeros(tTot2,1);
        RABPsave = zeros(tTot2,1);
    end
    
	RAinsave2 = zeros(tTot2,1);
	RAinsave3 = zeros(tTot2,1);
    RAinsave4 = zeros(tTot2,1);
	RARsave2  = zeros(tTot2,1);
	RARsave3  = zeros(tTot2,1);
	RARsave4  = zeros(tTot2,1);
    DiffX=DiffMat(alphaX,-2*alphaX,m,m,2,2,'x');
    DiffY=DiffMat(alphaY,-2*alphaY,n,n,2,2,'y');
    

	%GPU Setup
	if gpuEnabled
        X = gpuArray(X);
		RAout = gpuArray(RAout);
		RAin  = gpuArray(RAin);
		R     = gpuArray(R);
		RAR   = gpuArray(RAR);
		BP    = gpuArray(BP);
		RABP  = gpuArray(RABP);
		Hox   = gpuArray(Hox);
		Krox  = gpuArray(Krox);
		RAinsave = gpuArray(RAinsave);
		RARsave  = gpuArray(RARsave);
		production = gpuArray(production);
        RAoutNoise = gpuArray(RAoutNoise);
        RAinNoise  = gpuArray(RAinNoise);
        RARNoise   = gpuArray(RARNoise);
        HoxNoise   = gpuArray(HoxNoise);
        KroxNoise   = gpuArray(KroxNoise);
        
        RAoutTMP = gpuArray(RAoutTMP);
        RAinTMP  = gpuArray(RAinTMP);
        RTMP     = gpuArray(RTMP);
        RARTMP   = gpuArray(RARTMP);
        BPTMP    = gpuArray(BPTMP);
        RABPTMP  = gpuArray(RABPTMP);
        HoxTMP   = gpuArray(HoxTMP);
        KroxTMP  = gpuArray(KroxTMP);
        
        
		if initRun 
			RAoutsave= gpuArray(RAoutsave);
			Rsave    = gpuArray(Rsave);
			BPsave   = gpuArray(BPsave);
			RABPsave = gpuArray(RABPsave);
		end
		
		RAinsave2 = gpuArray(RAinsave2);
		RAinsave3 = gpuArray(RAinsave3);
		RAinsave4 = gpuArray(RAinsave4);
		RARsave2  = gpuArray(RARsave2);
		RARsave3  = gpuArray(RARsave3);
		RARsave4  = gpuArray(RARsave4);
        
        Vbp    = gpuArray(Vbp);
        gamma  = gpuArray(gamma);
        beta   = gpuArray(beta);
        kp     = gpuArray(kp);
        kmax   = gpuArray(kmax);
        beta0  = gpuArray(beta0);
        lambda = gpuArray(lambda);
        kdeg   = gpuArray(kdeg);
        vmax   = gpuArray(vmax);
        ron    = gpuArray(ron);
        roff   = gpuArray(roff);
        mon    = gpuArray(mon);
        moff   = gpuArray(moff);
        jalpha = gpuArray(jalpha);
        jbeta  = gpuArray(jbeta);
        bpdeg1 = gpuArray(bpdeg1);
        bpdeg2 = gpuArray(bpdeg2);
        Vr     = gpuArray(Vr);
        rdeg1  = gpuArray(rdeg1);
        rdeg2  = gpuArray(rdeg2);   
        d      = gpuArray(d);
        e      = gpuArray(e);
        cH     = gpuArray(cH);
        cK     = gpuArray(cK);
        kappaH = gpuArray(kappaH);
        kappaK = gpuArray(kappaK);
        dH     = gpuArray(dH);
        dK     = gpuArray(dK);
        DRA    = gpuArray(DRA);
        lambda = gpuArray(lambda);
        f0     = gpuArray(f0);
        
        epsinMult =gpuArray(epsinMult);
        epsinAdd  =gpuArray(epsinAdd);
        epsoutMult=gpuArray(epsoutMult);
        epsoutAdd =gpuArray(epsoutAdd);
        epsRAdd   =gpuArray(epsRAdd);
        epsRMult  =gpuArray(epsRMult);
        epsHAdd   =gpuArray(epsHAdd);
        epsHMult  =gpuArray(epsHMult);
        epsKAdd   =gpuArray(epsKAdd);
        epsKMult  =gpuArray(epsKMult);
		
		%Setup Diffusion Matrices
		DiffX=gpuArray(DiffX);
		DiffY=gpuArray(DiffY);
    else
        DiffX = sparse(DiffX);
        DiffY = sparse(DiffY);
    end

    maxState = 2 + enableState3;
    
    %Solver
    if initRun
        fprintf('\nEntering The Loop\n')
        progressbar
        loop = tic;
    end
    for state = 1:maxState
        if state == 2
            %Setup recording run
            if initRun
                progressbar
            end
            if bpMult==1 && cypMult==1
                gradMaxRAout = max(max(RAout));
                gradMaxRAin = max(max(RAin));
                gradMaxRAR = max(max(RAR));
                gradMinRAout = max(min(min(RAout)),0);
                gradMinRAin = max(min(min(RAin)),0);
                gradMinRAR = max(min(min(RAR)),0);
                gradSharpData = zeros(6,1);
                gradSharpData(1) = gradMaxRAout;
                gradSharpData(2) = gradMaxRAin;
                gradSharpData(3) = gradMaxRAR;
                gradSharpData(4) = gradMinRAout;
                gradSharpData(5) = gradMinRAin;
                gradSharpData(6) = gradMinRAR;
            else
                gradMaxRAout = gradSharpData(1);
                gradMaxRAin = gradSharpData(2);
                gradMaxRAR = gradSharpData(3);
                gradMinRAout = gradSharpData(4);
                gradMinRAin = gradSharpData(5);
                gradMinRAR = gradSharpData(6);
            end
            
            tTot = tTot2;
            dt = dt2;
            tEnd = tTot * dt;
            alpha  = sqrt(dt) / (dx*dx*dy*dy);
            
        end
        for i=1:tTot
            if initRun && mod(i,printSpacing)==0
                if benchmark
                    fprintf('%d\n',dt*i)
                end
                iter1 = tic;
                tic
            end
            %Solve Reactions
            if singlePerc && ~gpuEnabled
                DiffRA  = double(RAout)*DiffX + DiffY*double(RAout);
            else
                DiffRA  = RAout*DiffX + DiffY*RAout;
            end
			if gpuEnabled
					[RAoutN1,RAinN1,RN1,RARN1,BPN1,RABPN1,Hox1,Krox1] = arrayfun(@reactions,state,reaction,RAout,RAin,R,RAR,BP,RABP,Hox,Krox,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1, ...
                    Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK,RAoutTMP,RAinTMP,RTMP,RARTMP,BPTMP,RABPTMP,HoxTMP,KroxTMP);
					[RAoutN2,RAinN2,RN2,RARN2,BPN2,RABPN2,Hox2,Krox2] = arrayfun(@reactions,state,reaction,RAout + (dt/2)*RAoutN1,RAin+ (dt/2)*RAinN1,R + (dt/2)*RN1,RAR + (dt/2)*RARN1,BP+ (dt/2)*BPN1,RABP + (dt/2)*RABPN1,Hox + (dt/2)*Hox1,Krox + (dt/2)*Krox1,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK,RAoutTMP,RAinTMP,RTMP,RARTMP,BPTMP,RABPTMP,HoxTMP,KroxTMP);
					[RAoutN3,RAinN3,RN3,RARN3,BPN3,RABPN3,Hox3,Krox3] = arrayfun(@reactions,state,reaction,RAout + (dt/2)*RAoutN2,RAin+ (dt/2)*RAinN2,R + (dt/2)*RN2,RAR + (dt/2)*RARN2,BP+ (dt/2)*BPN2,RABP + (dt/2)*RABPN2,Hox + (dt/2)*Hox2,Krox + (dt/2)*Krox2,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK,RAoutTMP,RAinTMP,RTMP,RARTMP,BPTMP,RABPTMP,HoxTMP,KroxTMP);
					[RAoutN4,RAinN4,RN4,RARN4,BPN4,RABPN4,Hox4,Krox4] = arrayfun(@reactions,state,reaction,RAout + (dt)*RAoutN3,RAin+ (dt)*RAinN3,R + (dt)*RN3,RAR + (dt)*RARN3,BP+ (dt)*BPN3,RABP + (dt)*RABPN3,Hox + (dt)*Hox3,Krox + (dt)*Krox3,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK,RAoutTMP,RAinTMP,RTMP,RARTMP,BPTMP,RABPTMP,HoxTMP,KroxTMP);
            else
                    [RAoutN1,RAinN1,RN1,RARN1,BPN1,RABPN1,Hox1,Krox1] = reactions(state,reaction,RAout,RAin,R,RAR,BP,RABP,Hox,Krox,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK);
                %{
					[RAoutN2,RAinN2,RN2,RARN2,BPN2,RABPN2,Hox2,Krox2] = reactions(state,reaction,RAout + (dt/2)*RAoutN1,RAin+ (dt/2)*RAinN1,R + (dt/2)*RN1,RAR + (dt/2)*RARN1,BP+ (dt/2)*BPN1,RABP + (dt/2)*RABPN1,Hox + (dt/2)*Hox1,Krox + (dt/2)*Krox1,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK);
					[RAoutN3,RAinN3,RN3,RARN3,BPN3,RABPN3,Hox3,Krox3] = reactions(state,reaction,RAout + (dt/2)*RAoutN2,RAin+ (dt/2)*RAinN2,R + (dt/2)*RN2,RAR + (dt/2)*RARN2,BP+ (dt/2)*BPN2,RABP + (dt/2)*RABPN2,Hox + (dt/2)*Hox2,Krox + (dt/2)*Krox2,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK);
					[RAoutN4,RAinN4,RN4,RARN4,BPN4,RABPN4,Hox4,Krox4] = reactions(state,reaction,RAout + (dt)*RAoutN3,RAin+ (dt)*RAinN3,R + (dt)*RN3,RAR + (dt)*RARN3,BP+ (dt)*BPN3,RABP + (dt)*RABPN3,Hox + (dt)*Hox3,Krox + (dt)*Krox3,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1,...
					Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK);
                    %}
            end
            %Calculate Noise
            if state > 1
                if gpuEnabled
                    RAoutNoise = (epsoutAdd + epsoutMult*RAout ).*randn(n,m,'gpuArray')*alpha;
                    RAinNoise  = (epsinAdd  + epsinMult *RAin  ).*randn(n,m,'gpuArray')*alpha;
                    RARNoise   = (epsRAdd   + epsRMult  *RAR   ).*randn(n,m,'gpuArray')*alpha;
                    if state == 3
                        HoxNoise   = (epsHAdd   + epsHMult  *Hox   ).*randn(n,m,'gpuArray')*alpha;
                        KroxNoise  = (epsKAdd   + epsKMult  *Krox  ).*randn(n,m,'gpuArray')*alpha;
                    end
                else
                    RAoutNoise = (epsoutAdd + epsoutMult*RAout ).*randn(stream,n,m)*alpha;
                    RAinNoise  = (epsinAdd  + epsinMult *RAin  ).*randn(stream,n,m)*alpha;
                    RARNoise   = (epsRAdd   + epsRMult  *RAR   ).*randn(stream,n,m)*alpha;
                    if state == 3
                        HoxNoise   = (epsHAdd   + epsHMult  *Hox   ).*randn(stream,n,m)*alpha;
                        KroxNoise  = (epsKAdd   + epsKMult  *Krox  ).*randn(stream,n,m)*alpha;
                    end
                end
            end
            totChange = max(max(RAoutN1 + RAinN1 + RN1 + RARN1 + BPN1 + RABPN1));
            %Updates
            if gpuEnabled
                [RAout,RAin,R,RAR,BP,RABP,Hox,Krox] = arrayfun(@updateFn,state,dt,RAout,RAoutN1,RAoutN2,RAoutN3,RAoutN4,RAin,RAinN1,RAinN2,RAinN3,RAinN4,R,RN1,RN2,RN3,RN4,...
                RAR,RARN1,RARN2,RARN3,RARN4,BP,BPN1,BPN2,BPN3,BPN4,RABP,RABPN1,RABPN2,RABPN3,RABPN4,Hox,Hox1,Hox2,Hox3,Hox4,Krox,Krox1,Krox2,Krox3,Krox4,RAoutNoise,RAinNoise,RARNoise,HoxNoise,KroxNoise);
            else

                RAout = RAout + (dt)*(RAoutN1) + RAoutNoise ;
                RAin  = RAin  + (dt)*(RAinN1) + RAinNoise ;
                R     = R     + (dt)*(RN1);
                RAR   = RAR   + (dt)*(RARN1) + RARNoise ;
                BP    = BP    + (dt)*(BPN1);
                RABP  = RABP  + (dt)*(RABPN1);

               %{
                RAout = RAout + (dt/6)*(RAoutN1 + 2*RAoutN2 + 2*RAoutN3 + RAoutN4) + RAoutNoise ;
                RAin  = RAin  + (dt/6)*(RAinN1 + 2*RAinN2 + 2*RAinN3 + RAinN4) + RAinNoise ;
                R     = R     + (dt/6)*(RN1 + 2*RN2 + 2*RN3 + RN4);
                RAR   = RAR   + (dt/6)*(RARN1 + 2*RARN2 + 2*RARN3 + RARN4) + RARNoise ;
                BP    = BP    + (dt/6)*(BPN1 + 2*BPN2 + 2*BPN3 + BPN4);
                RABP  = RABP  + (dt/6)*(RABPN1 + 2*RABPN2 + 2*RABPN3 + RABPN4);
                %}
                if state == 3
                    Hox   = Hox   + (dt/6)*(Hox1 + 2*Hox2 + 2*Hox3 + Hox4)  + HoxNoise;
                    Krox  = Krox  + (dt/6)*(Krox1+ 2*Krox2+ 2*Krox3+ Krox4) + KroxNoise;
                end
            end
            if state == 2
                %Save
                RAinsave(i,1) =RAin(checkY,checkX1);
                RARsave(i,1)  = RAR(checkY,checkX1);
                RAinsave(i,2) =RAin(checkY,checkX2);
                RAinsave(i,3) =RAin(checkY,checkX3);
                RAinsave(i,4) =RAin(checkY,checkX4);
                RARsave(i,2)  = RAR(checkY,checkX2);
                RARsave(i,3)  = RAR(checkY,checkX3);
                RARsave(i,4)  = RAR(checkY,checkX4);

                if initRun
                    RAoutsave(i)=RAout(checkY,checkX4);
                    Rsave(i)    = R(checkY,checkX4);
                    BPsave(i)   = BP(checkY,checkX4);
                    RABPsave(i) = RABP(checkY,checkX4);
                end
            end
            if totChange < totChangeMin && state == 1 && i > 1000
                
                if initRun
                    fprintf('Total Change<%.2e => Steady State. Enter State 2\n',totChangeMin)
                end
                break;
            end
            if any(isnan(RAin(:)))>0
                RAin
                %NAN, Sim blew up. Throw error.
                error('SPDE:NAN','Error: RAin contains NAN. Sim Failed')
            end
            %Calculate Diffusion Terms        
            if initRun && mod(i,printSpacing)==0
                progressbar(i/tTot)
                if benchmark
                    fprintf('Reaction time: %s\n',toc(iter1))
                    %{
                    if state == 1
                        fprintf('Max RAout: %.2e\n',max(max(RAout)))
                        fprintf('Min RAout: %.2e\n',min(min(RAout)))
                    elseif state == 2
                        fprintf('Total Var RAin: %.2e\n',var(RAinsave4))
                        fprintf('Total Var RAR: %.2e\n',var(RARsave4))
                    end
                    %}
                    fprintf('Total Change: %.2e\n',totChange)
                end
            end
            if i==tTot && state == 1
                warning('A simulation prematurely exited State 1')
            end
        end
        if initRun
            toc(loop)
        end
    end
    
    %Calculate Sharpness
    for l = 1:length(sharpLocations)
        sigVal= (gradMaxRAout-gradMinRAout)*sharpLocations(l);
        sig50 = RAout > sigVal;
        sig50 = sig50 & cumsum(sig50,2) == 1;
        [Itemp,Jtemp] = find(sig50);
        Jtemp = Jtemp-ones(length(Jtemp),1);
        I = Itemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        J = Jtemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        if isempty(I) || isempty(J)
            BLRAout = nan;
        else
            %u0 = diag(RAout(I,J-ones(n,1)));
            u1 = diag(RAout(I,J));
            u2 = diag(RAout(I,J+ones(length(J),1)));
            %BLRAout = J*dx; %Zero Order
            interp = (sigVal*ones(length(J),1)-u2) ./ (u1 - u2);
            BLRAout = interp.*J.*dx + (ones(length(J),1)-interp).*(J+ones(length(J),1)).*dx; %Linear Interp Result
        end

        sigVal= (gradMaxRAin-gradMinRAin)*sharpLocations(l);
        sig50 = RAin > sigVal;
        sig50 = sig50 & cumsum(sig50,2) == 1;
        [Itemp,Jtemp] = find(sig50);
        Jtemp = Jtemp-ones(length(Jtemp),1);
        I = Itemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        J = Jtemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        if isempty(I) || isempty(J)
            BLRAin = nan;
        else
            %u0 = diag(RAin(I,J-ones(n,1)));
            u1 = diag(RAin(I,J));
            u2 = diag(RAin(I,J+ones(length(J),1)));
            %BLRAout = J*dx; %Zero Order
            interp = (sigVal*ones(length(J),1)-u2) ./ (u1 - u2);  
            BLRAin = interp.*J.*dx + (ones(length(J),1)-interp).*(J+ones(length(J),1)).*dx; %Linear Interp Result
        end

        sigVal= (gradMaxRAR-gradMinRAR)*sharpLocations(l);
        sig50 = RAR > sigVal;
        sig50 = sig50 & cumsum(sig50,2) == 1;
        [Itemp,Jtemp] = find(sig50);
        Jtemp = Jtemp-ones(length(Jtemp),1);
        I = Itemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        J = Jtemp(setdiff(1:length(Itemp),find(Jtemp==0)));
        if isempty(I) || isempty(J)
            BLRAR = nan;
        else
            %u0 = diag(RAR(I,J-ones(n,1)));
            u1 = diag(RAR(I,J));
            u2 = diag(RAR(I,J+ones(length(J),1)));
            %BLRAout = J*dx; %Zero Order
            interp = (sigVal*ones(length(J),1)-u2) ./ (u1 - u2);  
            BLRAR = interp.*J.*dx + (ones(length(J),1)-interp).*(J+ones(length(J),1)).*dx; %Linear Interp Result
        end

        sharpData((l-1)*6+1:l*6) = [mean(BLRAout),var(BLRAout),mean(BLRAin),var(BLRAin),mean(BLRAR),var(BLRAR)]; 
    end
    
	RAinGrad = RAin;
	RAoutGrad= RAout;
	RGrad    = R;
	BPGrad   = BP;
	RARGrad  = RAR;
	RABPGrad = RABP;
	HoxGrad  = Hox;
	KroxGrad = Krox;
	
	if gpuEnabled
		RAinGrad = gather(RAinGrad);
		RAoutGrad= gather(RAoutGrad);
		RGrad    = gather(RGrad);
		BPGrad   = gather(BPGrad);
		RARGrad  = gather(RARGrad);
		RABPGrad = gather(RABPGrad);
		HoxGrad  = gather(HoxGrad);
		KroxGrad = gather(KroxGrad);
		RAinsave = gather(RAinsave);
		RARsave  = gather(RARsave);
		
		if initRun 
			RAoutsave= gather(RAoutsave);
			Rsave    = gather(Rsave);
			BPsave   = gather(BPsave);
			RABPsave = gather(RABPsave);
		end
		
		RAinsave2 = gather(RAinsave2);
		RAinsave3 = gather(RAinsave3);
		RAinsave4 = gather(RAinsave4);
		RARsave2  = gather(RARsave2);
		RARsave3  = gather(RARsave3);
		RARsave4  = gather(RARsave4);
	end 

end