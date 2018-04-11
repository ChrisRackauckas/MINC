function [output] = driverSolo(description,dt,N,motif,varargin)
%Run Mode
id = [timeStamp() '-' motif];
name = 'RAsde';
dt     = str2num(dt);
N      = str2num(N);
motif  = str2num(motif);

parser = inputParser;
addOptional(parser,'gammaExpBias',0);
addOptional(parser,'deltaExpBias',0);
addOptional(parser,'lambdaExpBias',0);
addOptional(parser,'alphaExpBias',0);
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
alphaExpBias = parser.Results.alphaExpBias;

standardSetup(id,name,description,'doubleCores',batch);


betaBase   = 10^(-baseExp+2);
deltaBase  = 10^(-baseExp + deltaExpBias);
sigmaBase  = 10^(-baseExp+2);
sigma2Base = 10^(-baseExp+2);
sigma3Base = 10^(-baseExp+2);
etaBase    = 10^(-baseExp);
alpha0Base = 10^(-baseExp + alphaExpBias);
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

Fs = 1/dt;
f = 0:Fs/N:Fs/2;
figDefs = get(0,'defaultfigureposition');
diagnostics = false;
stream = RandStream('mrg32k3a');
tempRuns = true;

k=randi(50000);
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
    elseif motif ==23
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif23(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif3(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif31(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,RABP1] = simMotif32(vars(k,:),crabp,cypmax,N,dt,stream);
  elseif motif ==33
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,RABP1] = simMotif33(vars(k,:),crabp,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif4(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif41(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif42(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==43
		[RA,RAR,tmean1,tmean1R,tvar1,tvar1R,RABP1,R1,BP1] = simMotif43(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif5(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif51(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R] = simMotif52(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean1,tmean1R,tvar1,tvar1R,RABP1,R1,BP1] = simMotif53(vars(k,:),crabp,cypmax,N,dt,stream);
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
    
    RA1 = RA;
    RAR1 = RAR;
    
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
    elseif motif ==23
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif23(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif3(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif31(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,RABP2] = simMotif32(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif ==33
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,RABP2] = simMotif33(vars(k,:),crabpKnock,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif4(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif41(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif42(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==43
		[RA,RAR,tmean2,tmean2R,tvar2,tvar2R,RABP2,R2,BP2] = simMotif43(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif5(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif51(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R] = simMotif52(vars(k,:),crabpKnock,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean2,tmean2R,tvar2,tvar2R,RABP2,R2,BP2] = simMotif53(vars(k,:),crabpKnock,cypmax,N,dt,stream);
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
    
    RA2 = RA;
    RAR2 = RAR;
    
    cypmax = knock*cypmax;
    
    if motif ==1
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif1(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==2
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif2(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==21
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif21(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==22
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif22(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==23
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif23(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==3
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif3(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==31
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif31(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==32
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,RABP3] = simMotif32(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif ==33
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,RABP3] = simMotif33(vars(k,:),crabp,cypmax,N,dt,stream);
	elseif motif==4
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif4(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==41
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif41(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==42
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif42(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==43
		[RA,RAR,tmean3,tmean3R,tvar3,tvar3R,RABP3,R3,BP3] = simMotif43(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==5
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif5(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==51
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif51(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==52
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R] = simMotif52(vars(k,:),crabp,cypmax,N,dt,stream);
    elseif motif==53
        [RA,RAR,tmean3,tmean3R,tvar3,tvar3R,RABP3,R3,BP3] = simMotif53(vars(k,:),crabp,cypmax,N,dt,stream);
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
    
    RA3 = RA;
    RAR3 = RAR;
    
    perChangeMeanBP = (mean2-mean1)/max(mean1,mean2) * 100;
    perChangeVarBP  = (var2-var1)/max(var1,var2) * 100;
    perChangeCovBP  = (cov2-cov1)/max(cov1,cov2) * 100;
    
    perChangeMeanBPRAR = (mean2R-mean1R)/max(mean1R,mean2R) * 100;
    perChangeVarBPRAR  = (var2R-var1R)/max(var1R,var2R) * 100;
    perChangeCovBPRAR  = (cov2R-cov1R)/max(cov1R,cov2R) * 100;
    
    %Calculate mean and variance change
    perChangeMeanCyp    = (mean3-mean1)/max(mean1,mean3) * 100;
    perChangeVarCyp     = (var3-var1)/max(var1,var3) * 100;
    perChangeCovCyp     = (cov3-cov1)/max(cov1,cov3) * 100;
    perChangeMeanCypRAR = (mean3R-mean1R)/max(mean1R,mean3R) * 100;
    perChangeVarCypRAR  = (var3R-var1R)/max(var1R,var3R) * 100;
    perChangeCovCypRAR  = (cov3R-cov1R)/max(cov1R,cov3R) * 100;
    
    DeltaMeansBP=perChangeMeanBP;
    DeltaVarsBP =perChangeVarBP;
    DeltaCovsBP =perChangeCovBP;
    ZetasBP     =abs(perChangeVarBP)/(abs(perChangeMeanBP)+abs(perChangeVarBP));
    ZetasBP2    =abs(perChangeCovBP)/(abs(perChangeMeanBP)+abs(perChangeCovBP));
    
    DeltaMeansCyp=perChangeMeanCyp;
    DeltaVarsCyp =perChangeVarCyp;
    DeltaCovsCyp =perChangeCovCyp;
    ZetasCyp     =abs(perChangeVarCyp)/(abs(perChangeMeanCyp)+abs(perChangeVarCyp));
    ZetasCyp2    =abs(perChangeCovCyp)/(abs(perChangeMeanCyp)+abs(perChangeCovCyp));
    
    
   
    DeltaMeansBPRAR=perChangeMeanBPRAR;
    DeltaVarsBPRAR =perChangeVarBPRAR;
    DeltaCovsBPRAR =perChangeCovBPRAR;
    ZetasBPRAR     =abs(perChangeVarBPRAR)/(abs(perChangeMeanBPRAR)+abs(perChangeVarBPRAR));
    ZetasBPRAR2    =abs(perChangeCovBPRAR)/(abs(perChangeMeanBPRAR)+abs(perChangeCovBPRAR));
    
   
    DeltaMeansCypRAR =perChangeMeanCypRAR;
    DeltaVarsCypRAR  =perChangeVarCypRAR;
    DeltaCovsCypRAR  =perChangeCovCypRAR;
    ZetasCypRAR      =abs(perChangeVarCypRAR)/(abs(perChangeMeanCypRAR)+abs(perChangeVarCypRAR));
    ZetasCypRAR2     =abs(perChangeCovCypRAR)/(abs(perChangeMeanCypRAR)+abs(perChangeCovCypRAR));
    
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
varginTxt = true;
extendedScatter = false;
axisBot = -100;
edges = 0:.1:1;
%%
diaryTxt = '';
diaryTxt = standardFinish(id,name,description,dt,N,rndVar,motif,knock,batch);
if(batch)
    close(all)
end

%%

%{
disp(ZetasBP)
%disp(ZetasBP2)
disp(ZetasBPRAR)
%disp(ZetasBPRAR2)
disp(ZetasCyp)
%disp(ZetasCyp2)
disp(ZetasCypRAR)
%disp(ZetasCypRAR2)

tend = (N-1)*dt;
figure('Position',[figDefs(1),figDefs(2),1600,900])
subplot_tight(2,2,1,[margvScat marghScat]);
hold on
plot(0:dt:tend,RA1)
plot(0:dt:tend,RA2)
legend('Standard','BP Knockdown')
title(['RA ',num2str(ZetasBP)])
set(gca,'FontSize',fntSze)
hold off

subplot_tight(2,2,2,[margvScat marghScat]);
hold on
plot(0:dt:tend,RAR1)
plot(0:dt:tend,RAR2)
legend('Standard','BP Knockdown')
title(['RAR ',num2str(ZetasBPRAR)])
set(gca,'FontSize',fntSze)
hold off

subplot_tight(2,2,3,[margvScat marghScat]);
hold on
plot(0:dt:tend,RA1)
plot(0:dt:tend,RA3)
legend('Standard','Cyp Knockdown')
title(['RA ',num2str(ZetasCyp)])
set(gca,'FontSize',fntSze)
hold off

subplot_tight(2,2,4,[margvScat marghScat]);
hold on
plot(0:dt:tend,RAR1)
plot(0:dt:tend,RAR3)
legend('Standard','Cyp Knockdown')
title(['RAR ',num2str(ZetasCypRAR)])
set(gca,'FontSize',fntSze)
hold off

pubGCF('',10,'SDEKnockScatSSCov',varginTxt,id,email,diaryTxt,N,rndVar,motif)

%}
%{
tend = (N-1)*dt;
figure('Position',[figDefs(1),figDefs(2),1600,900])
subplot_tight(1,3,1,[margvScat marghScat]);
hold on
plot(0:dt:tend,RA1)
plot(0:dt:tend,RA2)
plot(0:dt:tend,RA3)

legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
title(['RA ',num2str(ZetasBP)])
set(gca,'FontSize',fntSze)
hold off

subplot_tight(1,3,2,[margvScat marghScat]);
hold on

plot(0:dt:tend,RABP1)
plot(0:dt:tend,RABP2)
plot(0:dt:tend,RABP3)

legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
title(['RAR ',num2str(ZetasBPRAR)])
set(gca,'FontSize',fntSze)
hold off

subplot_tight(1,3,3,[margvScat marghScat]);
hold on

plot(0:dt:tend,RAR1)
plot(0:dt:tend,RAR2)
plot(0:dt:tend,RAR3)
legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
title(['RA ',num2str(ZetasCyp)])
set(gca,'FontSize',fntSze)
hold off

pubGCF('',10,'SDEKnockScatSSCov',varginTxt,id,email,diaryTxt,N,rndVar,motif)
%}

if motif~=23
    tend = (N-1)*dt;
    figure('Position',[figDefs(1),figDefs(2),700,900])
    subplot_tight(3,1,1,[margvScat marghScat]);
    semilogy(0:dt:tend,RA1)
    hold on
    semilogy(0:dt:tend,RA2)
    semilogy(0:dt:tend,RA3)

    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['RA'])
    set(gca,'FontSize',fntSze)
    hold off

    subplot_tight(3,1,2,[margvScat marghScat]);
    semilogy(0:dt:tend,RABP1)
    hold on
    semilogy(0:dt:tend,RABP2)
    semilogy(0:dt:tend,RABP3)

    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['RABP'])
    set(gca,'FontSize',fntSze)
    hold off

    subplot_tight(3,1,3,[margvScat marghScat]);
    semilogy(0:dt:tend,RAR1)
    hold on
    semilogy(0:dt:tend,RAR2)
    semilogy(0:dt:tend,RAR3)
    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['RAR'])
    set(gca,'FontSize',fntSze)
    hold off

    pubGCF('',10,'TS',varginTxt,id,email,diaryTxt,N,rndVar,motif)
else
    
    tend = (N-1)*dt;
    figure('Position',[figDefs(1),figDefs(2),700,900])
    subplot_tight(2,1,1,[margvScat marghScat]);
    semilogy(0:dt:tend,RA1)
    hold on
    semilogy(0:dt:tend,RA2)
    semilogy(0:dt:tend,RA3)

    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['RA'])
    set(gca,'FontSize',fntSze)
    hold off

    subplot_tight(2,1,2,[margvScat marghScat]);
    semilogy(0:dt:tend,RAR1)
    hold on
    semilogy(0:dt:tend,RAR2)
    semilogy(0:dt:tend,RAR3)
    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['RAR'])
    set(gca,'FontSize',fntSze)
    hold off

    pubGCF('',10,'TS',varginTxt,id,email,diaryTxt,N,rndVar,motif)
end

if motif==43
    figure('Position',[figDefs(1),figDefs(2),700,600])
    subplot_tight(2,1,1,[margvScat marghScat]);
    semilogy(0:dt:tend,R1)
    hold on
    semilogy(0:dt:tend,R2)
    semilogy(0:dt:tend,R3)

    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['R'])
    set(gca,'FontSize',fntSze)
    hold off

    subplot_tight(2,1,2,[margvScat marghScat]);
    semilogy(0:dt:tend,BP1)
    hold on
    semilogy(0:dt:tend,BP2)
    semilogy(0:dt:tend,BP3)

    legend('Wildtype','BP-Knockdown','Cyp-Knockdown')
    title(['BP'])
    set(gca,'FontSize',fntSze)
    hold off

    pubGCF('',10,'TS2',varginTxt,id,email,diaryTxt,N,rndVar,motif)
end