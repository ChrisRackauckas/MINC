function [output] = driver(description,dt,N,K,motif,varargin)
%Run Mode
id = [timeStamp() '-' motif];
name = 'RAsde';
dt     = str2num(dt);
N      = str2num(N);
K      = str2num(K);
motif  = str2num(motif);

parser = inputParser;
addOptional(parser,'gammaExpBias',0);
addOptional(parser,'deltaExpBias',0);
addOptional(parser,'lambdaExpBias',0);
addOptional(parser,'nuExpBias',0);
addOptional(parser,'email',0);
addOptional(parser,'batch',0);
addOptional(parser,'singleComp',0);
addOptional(parser,'knock',.1);
addOptional(parser,'rndVar',5);
addOptional(parser,'baseExp',3);
parse(parser,varargin{:});
deltaExpBias = parser.Results.deltaExpBias;
email = parser.Results.email;
batch = parser.Results.batch;
singleComp = parser.Results.singleComp;
knock = parser.Results.knock;
rndVar = parser.Results.rndVar;
baseExp = parser.Results.baseExp;
gammaExpBias = parser.Results.gammaExpBias;
lambdaExpBias = parser.Results.lambdaExpBias;
nuExpBias = parser.Results.nuExpBias;

standardSetup(id,name,description,'doubleCores',batch);


betaBase   = 10^(-baseExp+1);
deltaBase  = 10^(-baseExp + deltaExpBias);
sigmaBase  = 10^(-baseExp);
sigma2Base = 10^(-baseExp);
sigma3Base = 10^(-baseExp);
etaBase    = 10^(-baseExp);
alpha0Base = 10^(-baseExp);
omega0Base = 10^(-baseExp);
gamma0Base = 10^(-baseExp + gammaExpBias);
lambda0Base= 10^(-baseExp + lambdaExpBias);
crabpBase  = 1;
cypmaxBase = 1;
nuBase     = 10^(-baseExp + nuExpBias);
rBase      = 10^(-baseExp);
zetaBase   = 10^(-baseExp);
aBase      = 10^(-baseExp);
uBase      = 10^(-baseExp);
bBase      = 10^(-baseExp);
cBase      = 10^(-baseExp);
dBase      = 10^(-baseExp);
eBase      = 10^(-baseExp);



calcFreq = false;

L = N/2 + 1;
DeltaMeansBP  = zeros(1,K);
DeltaVarsBP   = zeros(1,K);
ZetasBP       = zeros(1,K);
DeltaCovsBP   = zeros(1,K);
ZetasBP2      = zeros(1,K);

DeltaMeansCyp = zeros(1,K);
DeltaVarsCyp  = zeros(1,K);
ZetasCyp      = zeros(1,K);
DeltaCovsCyp  = zeros(1,K);
ZetasCyp2     = zeros(1,K);

DeltaMeansBPRAR  = zeros(1,K);
DeltaVarsBPRAR   = zeros(1,K);
ZetasBPRAR       = zeros(1,K);
DeltaCovsBPRAR   = zeros(1,K);
ZetasBPRAR2      = zeros(1,K);

DeltaMeansCypRAR = zeros(1,K);
DeltaVarsCypRAR  = zeros(1,K);
ZetasCypRAR      = zeros(1,K);
DeltaCovsCypRAR  = zeros(1,K);
ZetasCypRAR2     = zeros(1,K);

if motif == 7 || motif == 8 || motif == 9
    DeltaMeansBPD  = zeros(1,K);
    DeltaVarsBPD   = zeros(1,K);
    ZetasBPD       = zeros(1,K);
    DeltaCovsBPD   = zeros(1,K);
    ZetasBPD2       = zeros(1,K);
    
    DeltaMeansCypD = zeros(1,K);
    DeltaCovsCypD  = zeros(1,K);
    ZetasCypD2     = zeros(1,K);
end

meanBPDir        = zeros(1,K);
meanCypDir       = zeros(1,K);
varBPDir         = zeros(1,K);
varCypDir        = zeros(1,K);
covBPDir         = zeros(1,K);
covCypDir        = zeros(1,K);

vars             = zeros(K,20);
meanValues       = zeros(K,6);
varValues        = zeros(K,6);
covValues        = zeros(K,6);

if calcFreq
    freqSpect        = zeros(K,L);
    freqSpectBP      = zeros(K,L);
    freqSpectCyp     = zeros(K,L);
end
Fs = 1/dt;
f = 0:Fs/N:Fs/2;
figDefs = get(0,'defaultfigureposition');
diagnostics = false;
parfor_progress(id,K);
stream = RandStream('mrg32k3a');
tempRuns = false;

parfor k=1:K
    beta   = betaBase * randOoM(rndVar);
    delta  = deltaBase * randOoM(rndVar);
    sigma  = sigmaBase;
    sigma2 = sigma2Base;
    sigma3 = sigma3Base;
    cypmax = cypmaxBase * randOoM(rndVar);
    eta    = etaBase * randOoM(rndVar);
    alpha0 = alpha0Base * randOoM(rndVar);
    omega0 = omega0Base * randOoM(rndVar);
    gamma0 = gamma0Base * randOoM(rndVar);%-3);
    crabp  = crabpBase * randOoM(rndVar);

    nu     = nuBase * randOoM(rndVar);%-3);
    lambda0= lambda0Base * randOoM(rndVar);
    
    r      = rBase     * randOoM(rndVar);
    zeta   = zetaBase  * randOoM(rndVar);
    a      = aBase     * randOoM(rndVar);
    u      = uBase     * randOoM(rndVar);
	b      = bBase     * randOoM(rndVar);
    c      = cBase     * randOoM(rndVar);
    d      = dBase     * randOoM(rndVar);
    e      = eBase     * randOoM(rndVar);
    
    if singleComp
        vars(k,:) = single([beta,delta,sigma,cypmax,eta,alpha0,omega0,gamma0,nu,lambda0,r,zeta,a,u,b,c,d,e,sigma2,sigma3]);
    else
        vars(k,:) = [beta,delta,sigma,cypmax,eta,alpha0,omega0,gamma0,nu,lambda0,r,zeta,a,u,b,c,d,e,sigma2,sigma3];
    end
    
        
    set(stream,'Substream',k);
    
    if motif ==1
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif1(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==2
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif2(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==21
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif21(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==22
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif22(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif3(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif31(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif32(vars(k,:),crabp,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif4(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif41(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif42(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif5(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif51(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif52(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif53(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==54
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif54(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==55
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif55(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==56
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif56(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==57
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif57(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==58
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif58(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==6
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif6(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==7
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,D] = simMotif7(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==8
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,D] = simMotif8(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==9
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,D] = simMotif9(vars(k,:),crabp,cypmax,N,dt,stream);
    end
    
    mean1 = mean(RA);
    var1  = var(RA);
    mean1R= mean(RAR);
    var1R = var(RAR);
    cov1  = sqrt(var1)/mean1;
    cov1R = sqrt(var1R)/mean1R;
    
    tcov1 = sqrt(tvar1)/tmean1;
    tcov1R= sqrt(tvar1R)/tmean1R;
    
    if calcFreq
        xdft = fft(RA);
        xdft = xdft(1:N/2+1);
        P11 = (1/(Fs*N))*abs(xdft).^2;
        P11(2:end-1) = 2*P11(2:end-1);
    end
    
    if tempRuns
        RA1 = RA;
        RAR1 = RAR;
    end
    if motif == 7 || motif == 8 || motif == 9
        mean1D = mean(D);
        var1D  = var(D);
        cov1D  = sqrt(var1D)/mean1D;
    end
    
    crabpKnock= crabp*knock;
    
    %Reset RNG
    set(stream,'Substream',k);
    
    if motif ==1
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif1(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==2
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif2(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==21
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif21(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==22
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif22(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif3(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif31(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif32(vars(k,:),crabpKnock,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif4(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif41(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif42(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif5(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif51(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif52(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif53(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==54
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif54(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==55
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif55(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==56
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif56(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==57
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif57(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==58
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif58(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==6
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif6(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==7
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,D] = simMotif7(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==8
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,D] = simMotif8(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==9
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,D] = simMotif9(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    end
    
    tcov2 = sqrt(tvar2)/tmean2;
    tcov2R= sqrt(tvar2R)/tmean2R;
    
    mean2 = mean(RA);
    var2  = var(RA);
    mean2R= mean(RAR);
    var2R = var(RAR);
    cov2 = sqrt(var2)/mean2;
    cov2R= sqrt(var2R)/mean2R;
    
    if calcFreq
        xdft = fft(RA);
        xdft = xdft(1:N/2+1);
        P12 = (1/(Fs*N))*abs(xdft).^2;
        P12(2:end-1) = 2*P12(2:end-1);
    end
    
    if motif == 7 || motif == 8 || motif == 9
        mean2D = mean(D);
        var2D  = var(D);
        cov2D = sqrt(var2D)/mean2D;
    end
    if tempRuns
        RA2 = RA;
        RAR2 = RAR;
    end
    
    perChangeMeanBP = (mean2-mean1)/max(mean1,mean2) * 100;
    perChangeVarBP  = (var2-var1)/max(var1,var2) * 100;
    perChangeCovBP  = (cov2-cov1)/max(cov1,cov2) * 100;
    
    perChangeMeanBPRAR = (mean2R-mean1R)/max(mean1R,mean2R) * 100;
    perChangeVarBPRAR  = (var2R-var1R)/max(var1R,var2R) * 100;
    perChangeCovBPRAR  = (cov2R-cov1R)/max(cov1R,cov2R) * 100;
    
    if motif == 7 || motif == 8 || motif == 9
        perChangeMeanBPD = (mean2D-mean1D)/max(mean1D,mean2D) * 100;
        perChangeVarBPD  = (var2D-var1D)/max(var1D,var2D) * 100;
        perChangeCovBPD  = (cov2D-cov1D)/max(cov1D,cov2D) * 100;
    end
    
    perChangetMeanBP = (tmean2-tmean1)/max(tmean1,tmean2) * 100;
    perChangetVarBP  = (tvar2-tvar1)/max(tvar1,tvar2) * 100; 
    perChangetCovBP  = (tcov2-tcov1)/max(tcov1,tcov2) * 100; 
    perChangetMeanBPRAR = (tmean2R-tmean1R)/max(tmean1R,tmean2R) * 100;
    perChangetVarBPRAR = (tvar2R-tvar1R)/max(tvar1R,tvar2R) * 100;
    perChangetCovBPRAR  = (tcov2R-tcov1R)/max(tcov1R,tcov2R) * 100;
    
    cypmax= cypmax*knock;
    
    %Reset RNG
    set(stream,'Substream',k);
    
    if motif ==1
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif1(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==2
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif2(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==21
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif21(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==22
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif22(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif3(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif31(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif32(vars(k,:),crabp,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif4(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif41(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif42(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif5(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif51(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif52(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif53(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==54
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif54(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==55
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif55(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==56
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif56(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==57
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif57(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==58
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif58(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==6
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif6(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==7
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,D] = simMotif7(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==8
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,D] = simMotif8(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==9
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,D] = simMotif9(vars(k,:),crabp,cypmax,N,dt,stream);
    end
    
    tcov3 = sqrt(tvar3)/tmean3;
    tcov3R= sqrt(tvar3R)/tmean3R;
    
    mean3 = mean(RA);
    var3  = var(RA);
    mean3R= mean(RAR);
    var3R = var(RAR);
    
    cov3 = sqrt(var3)/mean3;
    cov3R = sqrt(var3R)/mean3R;
    
    if motif == 7 || motif == 8 || motif == 9
        mean3D = mean(D);
        var3D  = var(D);
        cov3D = sqrt(var3D)/mean3D;
    end
    
    if calcFreq
        xdft = fft(RA);
        xdft = xdft(1:N/2+1);
        P13 = (1/(Fs*N))*abs(xdft).^2;
        P13(2:end-1) = 2*P13(2:end-1);
    end
    
    if tempRuns
        RA3 = RA;
        RAR3 = RAR;
    end
    
    %Plot if not parfor
    if tempRuns
        tend = (N-1)*dt;
        figure
        hold on
        plot(0:dt:tend,RA1)
        plot(0:dt:tend,RA2)
        plot(0:dt:tend,RA3)
        legend('Standard','BP Knockdown','Cyp Knockdown')
        hold off

        if calcFreq
            figure
            loglog(f,smooth(P11,.001,'lowess'))
            hold on
            title('Single-Sided Amplitude Spectrum of X(t)')
            xlabel('f (Hz)')
            ylabel('|P1(f)|')
            loglog(f,smooth(P12,.001,'lowess'))
            loglog(f,smooth(P13,.001,'lowess'))
            legend('Standard','BP Knockdown','Cyp Knockdown')
            hold off

            figure
            freqDiff = abs(P13 - P11);
            loglog(f,smooth(freqDiff,.001,'lowess'))
            xlabel('f (Hz)')
            ylabel('|Wildtype(f) - CypKnockdown(f)|')

            figure
            freqDiff = abs(P12 - P11);
            loglog(f,smooth(freqDiff,.001,'lowess'))
            xlabel('f (Hz)')
            ylabel('|Wildtype(f) - BPKnockdown(f)|')
        end
    end
    
    %Calculate mean and variance change
    perChangeMeanCyp    = (mean3-mean1)/max(mean1,mean3) * 100;
    perChangeVarCyp     = (var3-var1)/max(var1,var3) * 100;
    perChangeCovCyp     = (cov3-cov1)/max(cov1,cov3) * 100;
    perChangeMeanCypRAR = (mean3R-mean1R)/max(mean1R,mean3R) * 100;
    perChangeVarCypRAR  = (var3R-var1R)/max(var1R,var3R) * 100;
    perChangeCovCypRAR  = (cov3R-cov1R)/max(cov1R,cov3R) * 100;
    
    if motif == 7 || motif == 8 || motif == 9 
        perChangeMeanCypD = (mean3D-mean1D)/max(mean1D,mean3D) * 100;
        perChangeVarCypD  = (var3D-var1D)/max(var1D,var3D) * 100;
        perChangeCovCypD  = (cov3D-cov1D)/max(cov1D,cov3D) * 100;
    end
    
    perChangetMeanCyp    = (tmean3-tmean1)/max(tmean1,tmean3) * 100;
    perChangetVarCyp     = (tvar3-tvar1)/max(tvar1,tvar3) * 100;
    perChangetCovCyp     = (tcov3-tcov1)/max(tcov1,tcov3) * 100;
    perChangetMeanCypRAR = (tmean3R-tmean1R)/max(tmean1R,tmean3R) * 100;
    perChangetVarCypRAR  = (tvar3R-tvar1R)/max(tvar1R,tvar3R) * 100;
    perChangetCovCypRAR  = (tcov3R-tcov1R)/max(tcov1R,tcov3R) * 100;
    
    DeltaMeansBP(k)=perChangeMeanBP;
    DeltaVarsBP(k) =perChangeVarBP;
    DeltaCovsBP(k) =perChangeCovBP;
    ZetasBP(k)     =abs(perChangeVarBP)/(abs(perChangeMeanBP)+abs(perChangeVarBP));
    ZetasBP2(k)    =abs(perChangeCovBP)/(abs(perChangeMeanBP)+abs(perChangeCovBP));
    
    DeltatMeansBP(k)=perChangetMeanBP;
    DeltatVarsBP(k) =perChangetVarBP;
    DeltatCovsBP(k) =perChangetCovBP;
    ZetastBP(k)     =abs(perChangetCovBP)/(abs(perChangetMeanBP)+abs(perChangetCovBP));
    ZetastBP2(k)    =abs(perChangetCovBP)/(abs(perChangetMeanBP)+abs(perChangetCovBP));
    
    DeltaMeansCyp(k)=perChangeMeanCyp;
    DeltaVarsCyp(k) =perChangeVarCyp;
    DeltaCovsCyp(k) =perChangeCovCyp;
    ZetasCyp(k)     =abs(perChangeVarCyp)/(abs(perChangeMeanCyp)+abs(perChangeVarCyp));
    ZetasCyp2(k)    =abs(perChangeCovCyp)/(abs(perChangeMeanCyp)+abs(perChangeCovCyp));
    
    DeltatMeansCyp(k)=perChangetMeanCyp;
    DeltatVarsCyp(k) =perChangetVarCyp;
    DeltatCovsCyp(k) =perChangetCovCyp;
    ZetastCyp(k)     =abs(perChangetVarCyp)/(abs(perChangetMeanCyp)+abs(perChangetVarCyp));
    ZetastCyp2(k)    =abs(perChangetCovCyp)/(abs(perChangetMeanCyp)+abs(perChangetCovCyp));
    
    DeltaMeansBPRAR(k)=perChangeMeanBPRAR;
    DeltaVarsBPRAR(k) =perChangeVarBPRAR;
    DeltaCovsBPRAR(k) =perChangeCovBPRAR;
    ZetasBPRAR(k)     =abs(perChangeVarBPRAR)/(abs(perChangeMeanBPRAR)+abs(perChangeVarBPRAR));
    ZetasBPRAR2(k)    =abs(perChangeCovBPRAR)/(abs(perChangeMeanBPRAR)+abs(perChangeCovBPRAR));
    
    DeltaMeansCypRAR(k)=perChangeMeanCypRAR;
    DeltaVarsCypRAR(k) =perChangeVarCypRAR;
    DeltaCovsCypRAR(k) =perChangeCovCypRAR;
    ZetasCypRAR(k)     =abs(perChangeVarCypRAR)/(abs(perChangeMeanCypRAR)+abs(perChangeVarCypRAR));
    ZetasCypRAR2(k)    =abs(perChangeCovCypRAR)/(abs(perChangeMeanCypRAR)+abs(perChangeCovCypRAR));
    
    DeltatMeansBPRAR(k)=perChangetMeanBPRAR;
    DeltatVarsBPRAR(k) =perChangetVarBPRAR;
    DeltatCovsBPRAR(k)  =perChangetCovBPRAR;
    ZetastBPRAR(k)     =abs(perChangetVarBPRAR)/(abs(perChangetMeanBPRAR)+abs(perChangetVarBPRAR));
    ZetastBPRAR2(k)    =abs(perChangetCovBPRAR)/(abs(perChangetMeanBPRAR)+abs(perChangetCovBPRAR));
    
    DeltatMeansCypRAR(k)=perChangetMeanCypRAR;
    DeltatVarsCypRAR(k) =perChangetVarCypRAR;
    DeltatCovsCypRAR(k) =perChangetCovCypRAR;
    ZetastCypRAR(k)     =abs(perChangetVarCypRAR)/(abs(perChangetMeanCypRAR)+abs(perChangetVarCypRAR));
    ZetastCypRAR2(k)    =abs(perChangetCovCypRAR)/(abs(perChangetMeanCypRAR)+abs(perChangetCovCypRAR));
    
    if motif == 7 || motif == 8 || motif == 9
        DeltaMeansBPD(k)=perChangeMeanBPD;
        DeltaVarsBPD(k) =perChangeVarBPD;
        DeltaCovsBPD(k) =perChangeCovBPD;
        ZetasBPD(k)     =abs(perChangeVarBPD)/(abs(perChangeMeanBPD)+abs(perChangeVarBPD));
        ZetasBPD2(k)    =abs(perChangeCovBPD)/(abs(perChangeMeanBPD)+abs(perChangeCovBPD));
        
        DeltaMeansCypD(k)=perChangeMeanCypD;
        DeltaVarsCypD(k) =perChangeVarCypD;
        DeltaCovsCypD(k) =perChangeCovCypD;
        ZetasCypD(k)     =abs(perChangeVarCypD)/(abs(perChangeMeanCypD)+abs(perChangeVarCypD));
        ZetasCypD2(k)    =abs(perChangeCovCypD)/(abs(perChangeMeanCypD)+abs(perChangeCovCypD));

    end
    
    %Store results
    
    meanBPDir(k)        = mean2>mean1;
    varBPDir(k)         = var2>var1;
    covBPDir(k)         = cov2>cov1;
    meanCypDir(k)       = mean3>mean1;
    varCypDir(k)        = var3>var1;
    covCypDir(k)        = cov3>cov1;
    
    meanBPDirR(k)        = mean2R>mean1R;
    varBPDirR(k)         = var2R>var1R;
    covBPDirR(k)         = cov2R>cov1R;
    meanCypDirR(k)       = mean3R>mean1R;
    varCypDirR(k)        = var3R>var1R;
    covCypDirR(k)        = cov3R>cov1R;
    
    if motif == 7 || motif == 8 || motif == 9
        meanBPDirD(k)        = mean2D>mean1D;
        varBPDirD(k)         = var2D>var1D;
        covBPDirD(k)         = cov2D>cov1D;
        meanCypDirD(k)       = mean3D>mean1D;
        varCypDirD(k)        = var3D>var1D;
        covCypDirR(k)        = cov3R>cov1R;
    end
    
    tmeanBPDir(k)        = tmean2>tmean1;
    tvarBPDir(k)         = tvar2>tvar1;
    tcovBPDir(k)         = tcov2>tcov1;
    tmeanCypDir(k)       = tmean3>tmean1;
    tvarCypDir(k)        = tvar3>tvar1;
    tcovCypDir(k)        = tcov3>tcov1;
    
    tmeanBPDirR(k)        = tmean2R>tmean1R;
    tvarBPDirR(k)         = tvar2R>tvar1R;
    tcovBPDirR(k)         = tcov2R>tcov1R;
    tmeanCypDirR(k)       = tmean3R>tmean1R;
    tvarCypDirR(k)        = tvar3R>tvar1R;
    tcovCypDirR(k)        = tcov3R>tcov1R;
    
    meanValues(k,:)     = [mean1 mean2 mean3 mean1R mean2R mean3R];
    varValues(k,:)      = [var1 var2 var3 var1R var2R var3R];
    covValues(k,:)      = [cov1 cov2 cov3 cov1R cov2R cov3R];
    
    tmeanValues(k,:)     = [tmean1 tmean2 tmean3 tmean1R tmean2R tmean3R];
    tvarValues(k,:)      = [tvar1 tvar2 tvar3 tvar1R tvar2R tvar3R];
    tcovValues(k,:)      = [tcov1 tcov2 tcov3 tcov1R tcov2R tcov3R];
    
    if calcFreq
        freqSpect(k,:)   = P11;
        freqSpectBP(k,:) = P12;
        freqSpectCyp(k,:)= P13;
    end
    parfor_progress(id);
end
parfor_progress(id,0);

%%
figDefs = get(0,'defaultfigureposition');
margvHist = .06;
marghHist = .1;
margvHistBot = .065;
padHist = .06;
margvScat = .07;
marghScat = .043;
padScat = .06;
fntSze = 15;
titSze = 15;
zetaSze = 25;
letterSize = 12;
histBins = 10;
map = [0 0 .5; .8 0 0];
varginTxt = false;
extendedScatter = false;
axisBot = -100;
edges = 0:.1:1;
%%
diaryTxt = '';
diaryTxt = standardFinish(id,name,description,dt,N,K,rndVar,motif,knock,batch);
if(batch)
    close(all)
end
delete(gcp)

axisBot = -100;
%% Mixed Plots
if motif ~= 7 && motif ~=8 && motif ~=9
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    %colormap(flipud(map))
    subplot_tight(2,1,1,[margvHist marghHist]);
    h1 = histogram(ZetasBP,edges,'FaceColor',[.8 0 0],'EdgeColor','black','FaceAlpha',1);
    hold on 
    h2 =histogram(ZetasCyp,edges,'FaceColor',[0 0 .5],'EdgeColor','black','FaceAlpha',.7);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Initial State',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(2,1,2,[margvHistBot marghHist]);
    h1 = histogram(ZetasBPRAR,edges);
    hold on 
    h2 =histogram(ZetasCypRAR,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title(sprintf('Motif %d: Terminal State',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
    set(gca,'FontSize',fntSze)
    hold off
    pubGCF('',10,'SDEKnockZeta',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    
    
    %Cov Plot
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(flipud(map))
    subplot_tight(2,1,1,[margvHist marghHist]);
    h1 = histogram(ZetasBP2,edges);
    hold on 
    h2 =histogram(ZetasCyp2,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Initial State',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(2,1,2,[margvHistBot marghHist]);
    h1 = histogram(ZetasBPRAR2,edges);
    hold on 
    h2 =histogram(ZetasCypRAR2,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title(sprintf('Motif %d: Terminal State',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
    set(gca,'FontSize',fntSze)
    hold off
    pubGCF('',10,'SDEKnockZeta2',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    %Scatter Plots Pooled
    figure('Position',[figDefs(1),figDefs(2),960,720])
    colormap(map)
    
    hold on
    subplot_tight(2,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltaMeansBP DeltaMeansCyp],[DeltaVarsBP DeltaVarsCyp],[ones(length(DeltaMeansBP),1)' zeros(length(DeltaMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Initial State, Cascade-Strength Knockdown',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(2,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltaMeansBPRAR DeltaMeansCypRAR],[DeltaVarsBPRAR DeltaVarsCypRAR],[ones(length(DeltaMeansBPRAR),1)' zeros(length(DeltaMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Terminal State, Cascade-Strength Knockdown',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    pubGCF('',10,'SDEKnockScat',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    %Scatter Plots CoV
    figure('Position',[figDefs(1),figDefs(2),960,720])
    colormap(map)
    
    hold on
    subplot_tight(2,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltaMeansBP DeltaMeansCyp],[DeltaCovsBP DeltaCovsCyp],[ones(length(DeltaMeansBP),1)' zeros(length(DeltaMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Initial State, Cascade-Strength Knockdown',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(2,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltaMeansBPRAR DeltaMeansCypRAR],[DeltaCovsBPRAR DeltaCovsCypRAR],[ones(length(DeltaMeansBPRAR),1)' zeros(length(DeltaMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Terminal State, Cascade-Strength Knockdown',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    pubGCF('',10,'SDEKnockScatCov',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    if motif < 4 || (motif > 20  && motif < 40)
        % Zeta Plots Theo
        
        figure('Position',[figDefs(1),figDefs(2),720,720])
    
        colormap(map)
        subplot_tight(2,1,1,[margvHist marghHist]);
        h1 = histogram(ZetastBP,edges);
        hold on 
        h2 =histogram(ZetastCyp,edges);
        h1.BinWidth = .10;
        h2.BinWidth = .10;
        axis([0 1 0 inf])
        title(sprintf('SS Motif %d: Initial State',motif))
        y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
        set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
        legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
        set(gca,'XTickLabel',[])
        set(gca,'FontSize',fntSze)
        hold off

        subplot_tight(2,1,2,[margvHist marghHist]);
        h1 = histogram(ZetastBPRAR,edges);
        hold on 
        h2 =histogram(ZetastCypRAR,edges);
        h1.BinWidth = .10;
        h2.BinWidth = .10;
        axis([0 1 0 inf])
        title(sprintf('SS Motif %d: Terminal State',motif))
        y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
        set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
        x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
        set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
        set(gca,'FontSize',fntSze)
        hold off

        pubGCF('',10,'SDEKnockZetaSS',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
        
        % Zeta Plots Cov
        
        figure('Position',[figDefs(1),figDefs(2),720,720])
    
        colormap(map)
        subplot_tight(2,1,1,[margvHist marghHist]);
        h1 = histogram(ZetastBP2,edges);
        hold on 
        h2 =histogram(ZetastCyp2,edges);
        h1.BinWidth = .10;
        h2.BinWidth = .10;
        axis([0 1 0 inf])
        title(sprintf('SS Motif %d: Initial State',motif))
        y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
        set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
        legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
        set(gca,'XTickLabel',[])
        set(gca,'FontSize',fntSze)
        hold off

        subplot_tight(2,1,2,[margvHist marghHist]);
        h1 = histogram(ZetastBPRAR2,edges);
        hold on 
        h2 =histogram(ZetastCypRAR2,edges);
        h1.BinWidth = .10;
        h2.BinWidth = .10;
        axis([0 1 0 inf])
        title(sprintf('SS Motif %d: Terminal State',motif))
        y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
        set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
        x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
        set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
        set(gca,'FontSize',fntSze)
        hold off

        pubGCF('',10,'SDEKnockZeta2SS',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
        
        
    %Scatter Plots SS
    figure('Position',[figDefs(1),figDefs(2),960,720])
    colormap(map)
    
    hold on
    subplot_tight(2,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltatMeansBP DeltatMeansCyp],[DeltatVarsBP DeltatVarsCyp],[ones(length(DeltatMeansBP),1)' zeros(length(DeltatMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Initial State SS',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(2,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltatMeansBPRAR DeltatMeansCypRAR],[DeltatVarsBPRAR DeltatVarsCypRAR],[ones(length(DeltatMeansBPRAR),1)' zeros(length(DeltatMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Terminal State SS',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    pubGCF('',10,'SDEKnockScatSS',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    %Scatter Plots CoV
    figure('Position',[figDefs(1),figDefs(2),960,720])
    colormap(map)
    
    hold on
    subplot_tight(2,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltatMeansBP DeltatMeansCyp],[DeltatCovsBP DeltatCovsCyp],[ones(length(DeltatMeansBP),1)' zeros(length(DeltatMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Initial State SS',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(2,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltatMeansBPRAR DeltatMeansCypRAR],[DeltatCovsBPRAR DeltatCovsCypRAR],[ones(length(DeltatMeansBPRAR),1)' zeros(length(DeltatMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Terminal State SS',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    pubGCF('',10,'SDEKnockScatSSCov',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    end
else %motif ==7 || 8
    
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(map)
    subplot_tight(3,1,1,[margvHist marghHist]);
    h1 = histogram(ZetasBP,edges);
    hold on 
    h2 =histogram(ZetasCyp,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Hematopoetic Stem Cells',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(3,1,2,[margvHist marghHist]);
    h1 = histogram(ZetasBPRAR,edges);
    hold on 
    h2 =histogram(ZetasCypRAR,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Lymphoid Primed',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(3,1,3,[margvHist marghHist]);
    h1 = histogram(ZetasBPD,edges);
    hold on 
    h2 =histogram(ZetasCypD,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title(sprintf('Motif %d: Myloid Primed',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
    set(gca,'FontSize',fntSze)
    hold off
    
    pubGCF('',10,'SDEKnockZeta',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    % CoV Zeta
   figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(map)
    subplot_tight(3,1,1,[margvHist marghHist]);
    h1 = histogram(ZetasBP2,edges);
    hold on 
    h2 =histogram(ZetasCyp2,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Hematopoetic Stem Cells',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(3,1,2,[margvHist marghHist]);
    h1 = histogram(ZetasBPRAR2,edges);
    hold on 
    h2 =histogram(ZetasCypRAR2,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    title(sprintf('Motif %d: Lymphoid Primed',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown','location','north');
    set(gca,'XTickLabel',[])
    set(gca,'FontSize',fntSze)
    hold off
    
    subplot_tight(3,1,3,[margvHist marghHist]);
    h1 = histogram(ZetasBPD2,edges);
    hold on 
    h2 =histogram(ZetasCypD2,edges);
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title(sprintf('Motif %d: Myloid Primed',motif))
    y = ylabel(sprintf('# of Simulations with \\fontsize{%d}\\zeta_{CoV}',fntSze));
    set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0]);
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',zetaSze),'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
    set(gca,'FontSize',fntSze)
    hold off
    
    pubGCF('',10,'SDEKnockZeta2',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    %Scatter Plots
    figure('Position',[figDefs(1),figDefs(2),960,720*1.5])
    colormap(map)
    hold on
    subplot_tight(3,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltaMeansBP DeltaMeansCyp],[DeltaVarsBP DeltaVarsCyp],[ones(length(DeltaMeansBP),1)' zeros(length(DeltaMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Hematopoetic Stem Cells',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(3,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltaMeansBPRAR DeltaMeansCypRAR],[DeltaVarsBPRAR DeltaVarsCypRAR],[ones(length(DeltaMeansBPRAR),1)' zeros(length(DeltaMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Lymphoid Primed',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(3,1,3,[margvScat marghScat]);
    colormap(map)
        gscatter([DeltaMeansBPD DeltaMeansCypD],[DeltaVarsBPD DeltaVarsCypD],[ones(length(DeltaMeansBPD),1)' zeros(length(DeltaMeansBPD),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Myloid Primed',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    pubGCF('',10,'SDEKnockScat',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    %Scatter Plots, Cov
    figure('Position',[figDefs(1),figDefs(2),960,720*1.5])
    margv = .07;
    margh = .043;
    pad = .06;
    colormap(map)
    hold on
    subplot_tight(3,1,1,[margvScat marghScat]);
    box on
    gscatter([DeltaMeansBP DeltaMeansCyp],[DeltaCovsBP DeltaCovsCyp],[ones(length(DeltaMeansBP),1)' zeros(length(DeltaMeansBP),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Hematopoetic Stem Cells',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(3,1,2,[margvScat marghScat]);
    colormap(map)
    gscatter([DeltaMeansBPRAR DeltaMeansCypRAR],[DeltaCovsBPRAR DeltaCovsCypRAR],[ones(length(DeltaMeansBPRAR),1)' zeros(length(DeltaMeansBPRAR),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Lymphoid Primed',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off
    
    hold on
    subplot_tight(3,1,3,[margvScat marghScat]);
    colormap(map)
        gscatter([DeltaMeansBPD DeltaMeansCypD],[DeltaCovsBPD DeltaCovsCypD],[ones(length(DeltaMeansBPD),1)' zeros(length(DeltaMeansBPD),1)'],map)
    axis([axisBot 100 axisBot 100])
    title(sprintf('Motif %d: Myloid Primed',motif))
    y = ylabel('%\Delta Variance','FontSize',fntSze);
    x = xlabel('%\Delta Mean','FontSize',fntSze);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown');
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',letterSize,'fontw','b')
    hold off    
    pubGCF('',10,'SDEKnockScatCov',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
end

%% Frequency Plots

if calcFreq
    figure('Position',[figDefs(1),figDefs(2),960,720*1.5])
    loglog(f,mean(freqSpect))
    hold on
    loglog(f,mean(freqSpectBP))
    loglog(f,mean(freqSpectCyp))
    title('Single-Sided Amplitude Spectrum of Initial State over Time')
    xlabel('f (Hz)')
    ylabel('|Init(f)|')
    legend('Standard','BP Knockdown','Cyp Knockdown')
    hold off
    pubGCF('',10,'meanFreqSpecs',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    figure('Position',[figDefs(1),figDefs(2),960,720*1.5])
    freqDiffBP = freqSpect - freqSpectBP;
    freqDiffCyp = freqSpect - freqSpectCyp;
    loglog(f,mean(freqDiffBP))
    xlabel('f (Hz)')
    ylabel('Wildtype(f) - BPKnockdown(f)')
    pubGCF('',10,'meanFreqSpecDiffBP',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
    
    figure('Position',[figDefs(1),figDefs(2),960,720*1.5])
    loglog(f,mean(freqDiffCyp))
    xlabel('f (Hz)')
    ylabel('|Wildtype(f) - CypKnockdown(f)|')
    pubGCF('',10,'meanFreqSpecDiffCyp',varginTxt,id,email,diaryTxt,N,K,rndVar,motif)
end

end

