function [RAtime,RARtime,RAinGrad,RARGrad,RAinGradBPknock,RARGradBPknock]=driver(description,tTot2,dt2,K,sigma,noiseRegime,vbpAdjust,enableState3,varargin)
id = timeStamp();
name = 'RAspde';

%Input Parse
parser = inputParser;

addOptional(parser,'tEnd',5e5);
addOptional(parser,'dt',5e-3);
addOptional(parser,'reaction',1);
addOptional(parser,'regime','Flat');
addOptional(parser,'diffOn',0);
addOptional(parser,'singlePerc',0);
addOptional(parser,'verbose',0);
addOptional(parser,'email',0);
addOptional(parser,'batch',0);
addOptional(parser,'gpus',0);
addOptional(parser,'benchmark',0);
addOptional(parser,'knockCyp',0);
addOptional(parser,'initRun',0);
addOptional(parser,'knock',.1);
addOptional(parser,'counterMax',8);
addOptional(parser,'errorDecMult',.8);
addOptional(parser,'dx',5);
addOptional(parser,'dy',5);
addOptional(parser,'totChangeMin',1e-10);
addOptional(parser,'firstRNGSeed',5);
addOptional(parser,'suppliedP',struct());
addOptional(parser,'useInits',true);

parse(parser,varargin{:});
tEnd = parser.Results.tEnd;
dt = parser.Results.dt;
reaction = parser.Results.reaction;
regime = parser.Results.regime;
diffOn = parser.Results.diffOn;
singlePerc = parser.Results.singlePerc;
verbose = parser.Results.verbose;
email = parser.Results.email;
batch = parser.Results.batch;
gpus = parser.Results.gpus;
benchmark = parser.Results.benchmark;
knock = parser.Results.knock;
counterMax = parser.Results.counterMax;
errorDecMult = parser.Results.errorDecMult;
dx = parser.Results.dx;
dy = parser.Results.dy;
totChangeMin = parser.Results.totChangeMin;
firstRNGSeed = parser.Results.firstRNGSeed;
knockCyp = parser.Results.knockCyp;
initRun = parser.Results.initRun;
suppliedP = parser.Results.suppliedP;
useInits = parser.Results.useInits;

warning('off','backtrace')

standardSetup(id,name,description);

xStart = -100;
xEnd   = 400;
yEnd   = 50;
m = (xEnd-xStart) / dx + 1; %Starts at 0
n = yEnd / dy +1 ; %Starts at 0
tTot = tEnd / dt + 1;
tEnd2 = tTot2*dt2;
sharpLocations = [.05 .1 .15 .2 .4 .6 .8]; %Manually change Blocation loop length, parfor error




%Run All Models
%Standard is taken at the middle of the gradient
%Taken at 1/5 along the x. i.e. 0, 100, 200, 300
%Non-RA all at 300

RAtime = zeros(tTot2,4,3); %dim1 = time, dim2 = location, dim3 = knock experiment (1 = none, 2 = bpknock, 3 = cypknock)
RARtime = zeros(tTot2,4,3);

gpuEnabled=gpus>=1;

if gpuEnabled
    gpuDevice(1);
    st = parallel.gpu.RandStream('Philox4x32-10','Seed',5);
    parallel.gpu.RandStream.setGlobalStream(st);
end
dtRun = dt;
preCounter = 0;

stream = RandStream('mlfg6331_64','Seed',firstRNGSeed);

if isempty(fieldnames(suppliedP)) %random parameters if not chosen
    p = parameterGen(sigma,regime,noiseRegime,diffOn,singlePerc,stream);
else
    p = suppliedP;
    if ~diffOn
        p.DRA = 0;
    end
end

if initRun
    %Initial Run
    reset(stream);
    while preCounter < counterMax
        try
            [RAtime(:,:,1),RARtime(:,:,1),sharpDataInit,gradSharpData,RAouttime,Rtime,BPtime,RABPtime,RAinGrad,RAoutGrad,RARGrad,RGrad,BPGrad,RABPGrad,HoxGrad,KroxGrad]      = main(tEnd,dx,dy,dtRun,1,1,p,reaction,true,tTot2,dt2,stream,gpuEnabled,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,nan,sharpLocations,useInits,knock);
            break;
        catch ME
            switch ME.identifier
                case 'SPDE:NAN'
                    dtRun = errorDecMult*dtRun;
                    preCounter = preCounter + 1;
                    warning('SPDE:DEC',sprintf('A decrement occured. dtRun=%e\n',dtRun))
                otherwise
                    rethrow(ME)
            end
        end
    end
    %BP Knockdown Run
    if gpuEnabled
        st = parallel.gpu.RandStream('Philox4x32-10','Seed',5);
        parallel.gpu.RandStream.setGlobalStream(st);
    end
    reset(stream);
    while preCounter < counterMax
        try
            [RAtime(:,:,2),RARtime(:,:,2),sharpDataBPKnockInit,~,RAoutBPknocktime,RBPknocktime,BPBPknocktime,RABPBPknocktime,RAinGradBPknock,RAoutGradBPknock,RARGradBPknock,RGradBPknock,BPGradBPknock,RABPGradBPknock,HoxGradBPknock,KroxGradBPknock] = main(tEnd,dx,dy,dtRun,1*knock,1,p,reaction,true,tTot2,dt2,stream,gpuEnabled,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,gradSharpData,sharpLocations,useInits,knock);
            break;
        catch ME
            switch ME.identifier
                case 'SPDE:NAN'
                    dtRun = errorDecMult*dtRun;
                    preCounter = preCounter + 1;
                    warning('SPDE:DEC',sprintf('A decrement occured. dtRun=%e\n',dtRun))
                otherwise
                    rethrow(ME)
            end
        end
    end

    if knockCyp
        if gpuEnabled
            st = parallel.gpu.RandStream('Philox4x32-10','Seed',5);
            parallel.gpu.RandStream.setGlobalStream(st);
        end
        reset(stream);
        while preCounter < counterMax
            try
                [RAtime(:,:,3),RARtime(:,:,3),sharpDataCypKnockInit,~,RAoutCypknocktime,RCypknocktime,BPCypknocktime,RABPCypknocktime,RAinGradCypknock,RAoutGradCypknock,RARGradCypknock,RGradCypknock,BPGradCypknock,RABPGradCypknock,HoxGradCypknock,KroxGradCypknock] = main(tEnd,dx,dy,dtRun,1,1*knock,p,reaction,true,tTot2,dt2,stream,gpuEnabled,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,gradSharpData,sharpLocations,useInits,knock);
                break;
            catch ME
                switch ME.identifier
                    case 'SPDE:NAN'
                        dtRun = errorDecMult*dtRun;
                        preCounter = preCounter + 1;
                        warning('A decrement occured. dtRun=%e\n',dtRun)
                    otherwise
                        rethrow(ME)
                end
            end
        end
    end

    if gpuEnabled
        gpuDevice([ ]);
    end

    fprintf('Diagnostics\n')
    fprintf('No Knockdown\n')
    fprintf('Mean and variance of RAout are: %.2e %.2e\n',mean(RAouttime),var(RAouttime))
    fprintf('Mean and variance of RAin are:  %.2e %.2e\n',mean(RAtime(:,4,1)),  var(RAtime(:,4,1)  ))
    fprintf('Mean and variance of RABP are:  %.2e %.2e\n',mean(RABPtime), var(RABPtime ))
    fprintf('Mean and variance of RAR are:   %.2e %.2e\n',mean(RARtime(:,4,1)), var(RARtime(:,4,1) ))
    fprintf('Mean and variance of R are:     %.2e %.2e\n',mean(Rtime),    var(Rtime    ))
    fprintf('Mean and variance of BP are:    %.2e %.2e\n',mean(BPtime),   var(BPtime   ))
    fprintf('CRABP Knockdown\n')
    fprintf('Mean and variance of RAout are: %.2e %.2e\n',mean(RAoutBPknocktime),var(RAoutBPknocktime))
    fprintf('Mean and variance of RAin are:  %.2e %.2e\n',mean(RAtime(:,4,2)),  var(RAtime(:,4,2)  ))
    fprintf('Mean and variance of RABP are:  %.2e %.2e\n',mean(RABPBPknocktime), var(RABPBPknocktime ))
    fprintf('Mean and variance of RAR are:   %.2e %.2e\n',mean(RARtime(:,4,2)), var(RARtime(:,4,2) ))
    fprintf('Mean and variance of R are:     %.2e %.2e\n',mean(RBPknocktime),    var(RBPknocktime    ))
    fprintf('Mean and variance of BP are:    %.2e %.2e\n',mean(BPBPknocktime),   var(BPBPknocktime   ))
    if knockCyp
        fprintf('Cyp Knockdown\n')
        fprintf('Mean and variance of RAout are: %.2e %.2e\n',mean(RAoutCypknocktime),var(RAoutCypknocktime))
        fprintf('Mean and variance of RAin are:  %.2e %.2e\n',mean(RAtime(:,4,3)),  var(RAtime(:,4,3)  ))
        fprintf('Mean and variance of RABP are:  %.2e %.2e\n',mean(RABPCypknocktime), var(RABPCypknocktime ))
        fprintf('Mean and variance of RAR are:   %.2e %.2e\n',mean(RARtime(:,4,3)), var(RARtime(:,4,3) ))
        fprintf('Mean and variance of R are:     %.2e %.2e\n',mean(RCypknocktime),    var(RCypknocktime    ))
        fprintf('Mean and variance of BP are:    %.2e %.2e\n',mean(BPCypknocktime),   var(BPCypknocktime   ))
    end
end

numFail = 0;
if K>0
    %Parallel Runs
    fprintf('\nStarting Parallel For Loop\n')
    stream = RandStream('mlfg6331_64');

    DeltaMeansBP  = zeros(4,K);
    DeltaVarsBP   = zeros(4,K);
    ZetasBP       = zeros(4,K);

    DeltaMeansCyp = zeros(4,K);
    DeltaVarsCyp  = zeros(4,K);
    ZetasCyp      = zeros(4,K);

    DeltaMeansBPRAR  = zeros(4,K);
    DeltaVarsBPRAR   = zeros(4,K);
    ZetasBPRAR       = zeros(4,K);

    DeltaMeansCypRAR = zeros(4,K);
    DeltaVarsCypRAR  = zeros(4,K);
    ZetasCypRAR      = zeros(4,K);

    meanBPDir        = zeros(4,K);
    meanCypDir       = zeros(4,K);
    varBPDir         = zeros(4,K);
    varCypDir        = zeros(4,K);

    meanBPDirR        = zeros(4,K);
    meanCypDirR       = zeros(4,K);
    varBPDirR         = zeros(4,K);
    varCypDirR        = zeros(4,K);

    parforTimer = tic;
    parfor_progressf(id,parforTimer,K);

    BDeltaMeansBP  = zeros(length(sharpLocations),K);
    BDeltaVarsBP   = zeros(length(sharpLocations),K);
    BZetasBP       = zeros(length(sharpLocations),K);

    BDeltaMeansCyp = zeros(length(sharpLocations),K);
    BDeltaVarsCyp  = zeros(length(sharpLocations),K);
    BZetasCyp      = zeros(length(sharpLocations),K);

    BDeltaMeansBPRAR  = zeros(length(sharpLocations),K);
    BDeltaVarsBPRAR   = zeros(length(sharpLocations),K);
    BZetasBPRAR       = zeros(length(sharpLocations),K);

    BDeltaMeansCypRAR = zeros(length(sharpLocations),K);
    BDeltaVarsCypRAR  = zeros(length(sharpLocations),K);
    BZetasCypRAR      = zeros(length(sharpLocations),K);

    BmeanBPDir        = zeros(length(sharpLocations),K);
    BmeanCypDir       = zeros(length(sharpLocations),K);
    BvarBPDir         = zeros(length(sharpLocations),K);
    BvarCypDir        = zeros(length(sharpLocations),K);

    BDeltaMeansBPX     = zeros(length(sharpLocations),K);
    BDeltaMeansCypX    = zeros(length(sharpLocations),K);
    BDeltaMeansBPXRAR  = zeros(length(sharpLocations),K);
    BDeltaMeansCypXRAR = zeros(length(sharpLocations),K);

    BmeanBPDirR        = zeros(length(sharpLocations),K);
    BmeanCypDirR       = zeros(length(sharpLocations),K);
    BvarBPDirR         = zeros(length(sharpLocations),K);
    BvarCypDirR        = zeros(length(sharpLocations),K);
    parforOff = false;

    parfor k=1:K
        warning('off','backtrace')
        gpuEnabled = gpus>=1;
        if gpuEnabled
            gpuID = selectGPUDevice();
            gpuEnabled = gpuID>=1;
        else
            gpuID = 0;
        end
        counter = 0;
        %Make Parameters
        if gpuEnabled
            st = parallel.gpu.RandStream('Philox4x32-10','Seed',k);
            parallel.gpu.RandStream.setGlobalStream(st);
        end
        set(stream,'Substream',k);
        p = parameterGen(sigma,regime,noiseRegime,diffOn,singlePerc,stream);
        dtParRun = dt;

        %Standard Run
        if gpuEnabled
            st = parallel.gpu.RandStream('Philox4x32-10','Seed',k);
            parallel.gpu.RandStream.setGlobalStream(st);
        end
        set(stream,'Substream',k);
        

        while counter < counterMax
            try
                [RA,RAR,sharpData,gradSharpDataPF]      = main(tEnd,dx,dy,dtParRun,1,1,p,reaction,parforOff,tTot2,dt2,stream,gpuID,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,nan,sharpLocations,useInits,knock);
                break;
            catch ME
                switch ME.identifier
                    case 'SPDE:NAN'
                        dtParRun = errorDecMult*dtParRun;
                        counter = counter + 1;
                        warning('SPDE:DEC',sprintf('A decrement occured on process id %2d. dtRun=%e',k,dtParRun))
                    otherwise
                        disp(k)
                        rethrow(ME)
                end
            end
        end

        %BP Knockdown Run
        if gpuEnabled
            st = parallel.gpu.RandStream('Philox4x32-10','Seed',k);
            parallel.gpu.RandStream.setGlobalStream(st);
        end
        set(stream,'Substream',k);
        while counter < counterMax
            try
                [RABPknock,RARBPknock,sharpDataBPKnock,~] = main(tEnd,dx,dy,dtParRun,1*knock,1,p,reaction,parforOff,tTot2,dt2,stream,gpuID,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,gradSharpDataPF,sharpLocations,useInits,knock);
                break;
            catch ME
                switch ME.identifier
                    case 'SPDE:NAN'
                        dtParRun = errorDecMult*dtParRun;
                        counter = counter + 1;
                        warning('SPDE:DEC',sprintf('A decrement occured on process id %2d. dtRun=%e',k,dtParRun))
                    otherwise
                        disp(k)
                        rethrow(ME)
                end
            end
        end

        %Cyp Knockdown Run
        if knockCyp
            if gpuEnabled
                st = parallel.gpu.RandStream('Philox4x32-10','Seed',k);
                parallel.gpu.RandStream.setGlobalStream(st);
            end
            set(stream,'Substream',k);
            while counter < counterMax
                try
                    [RACypknock,RARCypknock,sharpDataCypKnock,~] = main(tEnd,dx,dy,dtParRun,1,1*knock,p,reaction,parforOff,tTot2,dt2,stream,gpuID,totChangeMin,singlePerc,vbpAdjust,enableState3,benchmark,gradSharpDataPF,sharpLocations,useInits,knock);
                    break;
                catch ME
                    switch ME.identifier
                        case 'SPDE:NAN'
                            dtParRun = errorDecMult*dtParRun;
                            counter = counter + 1;
                            warning('SPDE:DEC',sprintf('A decrement occured on process id %2d. dtRun=%e',k,dtParRun))
                        otherwise
                            disp(k)
                            rethrow(ME)
                    end
                end
            end
        end

        if gpuEnabled
            gpuDevice([ ]);
        end
        if counter >= counterMax
            warning('Run %d overdecremented and is excluded from results.',k)
            mean1 = nan(1,4);
            var1  = nan(1,4);
            mean1R= nan(1,4);
            var1R = nan(1,4);

            mean2 = nan(1,4);
            var2  = nan(1,4);
            mean2R= nan(1,4);
            var2R = nan(1,4);

            mean3 = nan(1,4);
            var3  = nan(1,4);
            mean3R= nan(1,4);
            var3R = nan(1,4);
            numFail = numFail + 1;

            Bmean1 = nan(1,length(sharpLocations));
            Bvar1  = nan(1,length(sharpLocations));
            Bmean1R= nan(1,length(sharpLocations));
            Bvar1R = nan(1,length(sharpLocations));

            Bmean2 = nan(1,length(sharpLocations));
            Bvar2  = nan(1,length(sharpLocations));
            Bmean2R= nan(1,length(sharpLocations));
            Bvar2R = nan(1,length(sharpLocations));

            Bmean3 = nan(1,length(sharpLocations));
            Bvar3  = nan(1,length(sharpLocations));
            Bmean3R= nan(1,length(sharpLocations));
            Bvar3R = nan(1,length(sharpLocations));
        else
            mean1 = mean(RA);
            var1  = var(RA);
            mean1R= mean(RAR);
            var1R = var(RAR);


            mean2 = mean(RABPknock);
            var2  = var(RABPknock);
            mean2R= mean(RARBPknock);
            var2R = var(RARBPknock);
            if knockCyp
                mean3 = mean(RACypknock);
                var3  = var(RACypknock);
                mean3R= mean(RARCypknock);
                var3R = var(RARCypknock);
            end

            Bmean1 = sharpData(3:6:6*length(sharpLocations));
            Bvar1  = sharpData(4:6:6*length(sharpLocations));
            Bmean1R= sharpData(5:6:6*length(sharpLocations));
            Bvar1R=  sharpData(6:6:6*length(sharpLocations));

            Bmean2 = sharpDataBPKnock(3:6:6*length(sharpLocations));
            Bvar2  = sharpDataBPKnock(4:6:6*length(sharpLocations));
            Bmean2R= sharpDataBPKnock(5:6:6*length(sharpLocations));
            Bvar2R=  sharpDataBPKnock(6:6:6*length(sharpLocations));
            if knockCyp
                Bmean3 = sharpDataCypKnock(3:6:6*length(sharpLocations));
                Bvar3  = sharpDataCypKnock(4:6:6*length(sharpLocations));
                Bmean3R= sharpDataCypKnock(5:6:6*length(sharpLocations));
                Bvar3R=  sharpDataCypKnock(6:6:6*length(sharpLocations));
            end
        end
        for location = 1:4
            perChangeMeanBP = abs((mean2(location)-mean1(location))/max(mean1(location),mean2(location)) * 100);
            perChangeVarBP  = abs((var2(location)-var1(location))/max(var1(location),var2(location)) * 100);
            perChangeMeanBPRAR = abs((mean2R(location)-mean1R(location))/max(mean1R(location),mean2R(location)) * 100);
            perChangeVarBPRAR  = abs((var2R(location)-var1R(location))/max(var1R(location),var2R(location)) * 100);

            DeltaMeansBP(location,k)=perChangeMeanBP;
            DeltaVarsBP(location,k) =perChangeVarBP;
            ZetasBP(location,k)     =perChangeVarBP/(perChangeMeanBP+perChangeVarBP);

            DeltaMeansBPRAR(location,k)=perChangeMeanBPRAR;
            DeltaVarsBPRAR(location,k) =perChangeVarBPRAR;
            ZetasBPRAR(location,k)     =perChangeVarBPRAR/(perChangeMeanBPRAR+perChangeVarBPRAR);

            meanBPDir(location,k)        = mean2(location)>mean1(location);
            varBPDir(location,k)         = var2(location)>var1(location);


            meanBPDirR(location,k)        = mean2R(location)>mean1R(location);
            varBPDirR(location,k)         = var2R(location)>var1R(location);

            if knockCyp
                perChangeMeanCyp = abs((mean3(location)-mean1(location))/max(mean1(location),mean3(location)) * 100);
                perChangeVarCyp  = abs((var3(location)-var1(location))/max(var1(location),var3(location)) * 100);
                perChangeMeanCypRAR = abs((mean3R(location)-mean1R(location))/max(mean1R(location),mean3R(location)) * 100);
                perChangeVarCypRAR  = abs((var3R(location)-var1R(location))/max(var1R(location),var3R(location)) * 100);
                
                DeltaMeansCyp(location,k)=perChangeMeanCyp;
                DeltaVarsCyp(location,k) =perChangeVarCyp;
                ZetasCyp(location,k)     =perChangeVarCyp/(perChangeMeanCyp+perChangeVarCyp);

                DeltaMeansCypRAR(location,k)=perChangeMeanCypRAR;
                DeltaVarsCypRAR(location,k) =perChangeVarCypRAR;
                ZetasCypRAR(location,k)     =perChangeVarCypRAR/(perChangeMeanCypRAR+perChangeVarCypRAR);
                
                meanCypDir(location,k)       = mean3(location)>mean1(location);
                varCypDir(location,k)        = var3(location)>var1(location);
                meanCypDirR(location,k)       = mean3R(location)>mean1R(location);
                varCypDirR(location,k)        = var3R(location)>var1R(location);
            end
        end
        for Blocation=1:7
            % Boundary Terms

            BperChangeMeanBP = abs((Bmean2(Blocation)-Bmean1(Blocation))/max(Bmean1(Blocation),Bmean2(Blocation)) * 100);
            BChangeMeanBP = abs((Bmean2(Blocation)-Bmean1(Blocation)));
            BperChangeVarBP  = abs((Bvar2(Blocation)-Bvar1(Blocation))/max(Bvar1(Blocation),Bvar2(Blocation)) * 100);
            BperChangeMeanBPRAR = abs((Bmean2R(Blocation)-Bmean1R(Blocation))/max(Bmean1R(Blocation),Bmean2R(Blocation)) * 100);
            BChangeMeanBPRAR = abs((Bmean2R(Blocation)-Bmean1R(Blocation)));
            BperChangeVarBPRAR  = abs((Bvar2R(Blocation)-Bvar1R(Blocation))/max(Bvar1R(Blocation),Bvar2R(Blocation)) * 100);

            BDeltaMeansBP(Blocation,k)=BperChangeMeanBP;
            BDeltaVarsBP(Blocation,k) =BperChangeVarBP;
            BZetasBP(Blocation,k)     =BperChangeVarBP/(BperChangeMeanBP+BperChangeVarBP);

            

            BDeltaMeansBPRAR(Blocation,k)=BperChangeMeanBPRAR;
            BDeltaVarsBPRAR(Blocation,k) =BperChangeVarBPRAR;
            BZetasBPRAR(Blocation,k)     =BperChangeVarBPRAR/(BperChangeMeanBPRAR+BperChangeVarBPRAR);



            BDeltaMeansBPX(Blocation,k)  = BChangeMeanBP;
            
            BDeltaMeansBPXRAR(Blocation,k)  = BChangeMeanBPRAR;
            

            BmeanBPDir(Blocation,k)        = Bmean2(Blocation)>Bmean1(Blocation);
            BvarBPDir(Blocation,k)         = Bvar2(Blocation)>Bvar1(Blocation);


            BmeanBPDirR(Blocation,k)        = Bmean2R(Blocation)>Bmean1R(Blocation);
            BvarBPDirR(Blocation,k)         = Bvar2R(Blocation)>Bvar1R(Blocation);
            
            if knockCyp
                BperChangeMeanCyp = abs((Bmean3(Blocation)-Bmean1(Blocation))/max(Bmean1(Blocation),Bmean3(Blocation)) * 100);
                BChangeMeanCyp = abs((Bmean3(Blocation)-Bmean1(Blocation)));
                BperChangeVarCyp  = abs((Bvar3(Blocation)-Bvar1(Blocation))/max(Bvar1(Blocation),Bvar3(Blocation)) * 100);
                BperChangeMeanCypRAR = abs((Bmean3R(Blocation)-Bmean1R(Blocation))/max(Bmean1R(Blocation),Bmean3R(Blocation)) * 100);
                BChangeMeanCypRAR = abs((Bmean3R(Blocation)-Bmean1R(Blocation)));
                BperChangeVarCypRAR  = abs((Bvar3R(Blocation)-Bvar1R(Blocation))/max(Bvar1R(Blocation),Bvar3R(Blocation)) * 100);
                
                BDeltaMeansCyp(Blocation,k)=BperChangeMeanCyp;
                BDeltaVarsCyp(Blocation,k) =BperChangeVarCyp;
                BZetasCyp(Blocation,k)     =BperChangeVarCyp/(BperChangeMeanCyp+BperChangeVarCyp);
                BDeltaMeansCypRAR(Blocation,k)=BperChangeMeanCypRAR;
                BDeltaVarsCypRAR(Blocation,k) =BperChangeVarCypRAR;
                BZetasCypRAR(Blocation,k)     =BperChangeVarCypRAR/(BperChangeMeanCypRAR+BperChangeVarCypRAR);
                BDeltaMeansCypX(Blocation,k)  = BChangeMeanCyp;
                BDeltaMeansCypXRAR(Blocation,k)  = BChangeMeanCypRAR;
                BmeanCypDir(Blocation,k)       = Bmean3(Blocation)>Bmean1(Blocation);
                BvarCypDir(Blocation,k)        = Bvar3(Blocation)>Bvar1(Blocation);
                BmeanCypDirR(Blocation,k)       = Bmean3R(Blocation)>Bmean1R(Blocation);
                BvarCypDirR(Blocation,k)        = Bvar3R(Blocation)>Bvar1R(Blocation);
            end
        end

        parfor_progressf(id,parforTimer);
    end
    parfor_progressf(id,parforTimer,0);
    toc(parforTimer);
end

fprintf('\nLoop Ended.\n\nThe number of runs that failed is %d.\n',numFail)
%Set colormap
map = [0 0 .5; .8 0 0];

diaryTxt = '';
diaryTxt = standardFinish(id,name,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,firstRNGSeed,enableState3,verbose,email,batch,gpus);

if verbose
    if initRun
        if knockCyp
            iMax = 3;
        else
            iMax = 2;
        end
        %% Time Series Plots
        figure('Position',[100,100,960,720])
        subplot(2,3,1)
        hold on
        for i = 1:iMax
            plot(0:dt2:tEnd2-dt2,RAtime(:,4,i))
        end
        ylabel('RAin','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        subplot(2,3,2)
        hold on
        plot(0:dt2:tEnd2-dt2,RAouttime)
        plot(0:dt2:tEnd2-dt2,RAoutBPknocktime)
        if knockCyp
            plot(0:dt2:tEnd2-dt2,RAoutCypknocktime)
        end
        ylabel('RAout','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        subplot(2,3,3)
        hold on
        for i = 1:iMax
            plot(0:dt2:tEnd2-dt2,RARtime(:,4,i))
        end
        ylabel('RAR','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        subplot(2,3,4)
        hold on
        plot(0:dt2:tEnd2-dt2,RABPtime)
        plot(0:dt2:tEnd2-dt2,RABPBPknocktime)
        if knockCyp
            plot(0:dt2:tEnd2-dt2,RABPCypknocktime)
        end
        ylabel('RABP','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        subplot(2,3,5)
        hold on
        plot(0:dt2:tEnd2-dt2,Rtime)
        plot(0:dt2:tEnd2-dt2,RBPknocktime)
        if knockCyp
            plot(0:dt2:tEnd2-dt2,RCypknocktime)
        end
        ylabel('R','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        subplot(2,3,6)
        hold on
        plot(0:dt2:tEnd2-dt2,BPtime)
        plot(0:dt2:tEnd2-dt2,BPBPknocktime)
        if knockCyp
            plot(0:dt2:tEnd2-dt2,BPCypknocktime)
        end
        ylabel('BP','FontSize',9)
        xlabel('t','FontSize',9)
        hold off

        leg = legend('Standard','BP Knock','Cyp Knock');
        % Programatically move the Legend
        newPosition = [0.9 0.4 0.01 0.1];
        newUnits = 'normalized';
        set(leg,'Position', newPosition,'Units', newUnits,'FontSize',8);

        noteGCF('X=300 Time Series Plots',10,'SPDETS_Mid', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

        %% RA/RAR Time Series Plots
        figure('Position',[100,100,1280,720])
        for loc = 1:4
            subplot(2,4,loc)
            hold on
            for i = 1:iMax
                plot(0:dt2:tEnd2-dt2,RAtime(:,loc,i))
            end
            ylabel(sprintf('RAin%d',loc),'FontSize',9)
            xlabel('t','FontSize',9)
            hold off

            subplot(2,4,loc+4)
            hold on
            for i = 1:iMax
                plot(0:dt2:tEnd2-dt2,RARtime(:,loc,i))
            end
            ylabel(sprintf('RAR%d',loc),'FontSize',9)
            xlabel('t','FontSize',9)
            hold off
        end

        leg = legend('Standard','BP Knock','Cyp Knock');
        % Programatically move the Legend
        newPosition = [0.9 0.4 0.01 0.1];
        newUnits = 'normalized';
        set(leg,'Position', newPosition,'Units', newUnits,'FontSize',8);

        noteGCF('RA/RAR Time Series Plots',10,'SPDETS_Mid', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

        %% Gradient Plots
        checkY  = round(n/2); %Plot the gradient along x from the middle y
        figure('Position',[100,100,960,720])
        subplot(2,3,1)
        hold on
        plot(xStart:dx:xEnd,RAinGrad(checkY,:))
        plot(xStart:dx:xEnd,RAinGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,RAinGradCypknock(checkY,:))
        end
        ylabel('RAin','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        subplot(2,3,2)
        hold on
        plot(xStart:dx:xEnd,RAoutGrad(checkY,:))
        plot(xStart:dx:xEnd,RAoutGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,RAoutGradCypknock(checkY,:))
        end
        ylabel('RAout','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        subplot(2,3,3)
        hold on
        plot(xStart:dx:xEnd,RARGrad(checkY,:))
        plot(xStart:dx:xEnd,RARGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,RARGradCypknock(checkY,:))
        end
        ylabel('RAR','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        subplot(2,3,4)
        hold on
        plot(xStart:dx:xEnd,RABPGrad(checkY,:))
        plot(xStart:dx:xEnd,RABPGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,RABPGradCypknock(checkY,:))
        end
        ylabel('RABP','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        subplot(2,3,5)
        hold on
        plot(xStart:dx:xEnd,RGrad(checkY,:))
        plot(xStart:dx:xEnd,RGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,RGradCypknock(checkY,:))
        end
        ylabel('R','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        subplot(2,3,6)
        hold on
        plot(xStart:dx:xEnd,BPGrad(checkY,:))
        plot(xStart:dx:xEnd,BPGradBPknock(checkY,:))
        if knockCyp
            plot(xStart:dx:xEnd,BPGradCypknock(checkY,:))
        end
        ylabel('BP','FontSize',9)
        xlabel('x','FontSize',9)
        hold off

        leg = legend('Standard','BP Knock','Cyp Knock');
        % Programatically move the Legend
        newPosition = [0.9 0.4 0.01 0.1];
        newUnits = 'normalized';
        set(leg,'Position', newPosition,'Units', newUnits,'FontSize',8);


        noteGCF('End Time Gradients',20,'SPDEGrad', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

     
        %% 2D Gradient Plots
        XimgStart = 210;
        XimgEnd = 300;

        StartIdx = (XimgStart-xStart)/dx;
        EndIdx = length(xStart:dx:xEnd)-(xEnd-XimgEnd)/dx;
        
        figure('Position',[100,100,1480,720])
        subplot(2,3,1)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAinGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAin')
        xlabel('x','Fontsize',1)
        ylabel('y','Fontsize',1)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,2)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAoutGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAout')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,3)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RARGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAR')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,4)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RABPGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RABP')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,5)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, R')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,6)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,BPGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, BP')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        noteGCF('End Time Gradients',20,'SPDE2DGrad', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)
        

        StartIdx = (XimgStart-xStart)/dx;
        EndIdx = length(xStart:dx:xEnd)-(xEnd-XimgEnd)/dx;
        
        figure('Position',[100,100,1480,720])
        subplot(2,3,1)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAinGradBPknock(:,StartIdx:EndIdx,:))
        colorbar
        title('Gradient Plot, RAin')
        xlabel('x','Fontsize',1)
        ylabel('y','Fontsize',1)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,2)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAoutGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAout')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,3)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RARGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAR')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,4)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RABPGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RABP')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,5)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, R')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        subplot(2,3,6)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,BPGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, BP')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        hold off

        noteGCF('End Time Gradients',20,'SPDE2DGrad', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)
    
    %{
        %% CumVar Plots
        figure('Position',[100,100,1280,720])
        for loc = 1:4
            [m,v]   = cumMeanVar(RAtime(:,loc,1));
            [mBP,vBP]   = cumMeanVar(RAtime(:,loc,2));
            [mC,vC]   = cumMeanVar(RAtime(:,loc,3));
            [mR,vR]   = cumMeanVar(RARtime(:,loc,1));
            [mBPR,vBPR]   = cumMeanVar(RARtime(:,loc,2));
            [mCR,vCR]   = cumMeanVar(RARtime(:,loc,3));

            subplot(2,4,loc)
            hold on
            plot(0:dt2:tEnd2-dt2,v)
            plot(0:dt2:tEnd2-dt2,vBP)
            plot(0:dt2:tEnd2-dt2,vC)
            ylabel(sprintf('RAin%d Cumulative Variance',loc),'FontSize',9)
            xlabel('t','FontSize',9)
            hold off

            subplot(2,4,loc+4)
            hold on
            plot(0:dt2:tEnd2-dt2,vR)
            plot(0:dt2:tEnd2-dt2,vBPR)
            plot(0:dt2:tEnd2-dt2,vCR)
            ylabel(sprintf('RAR%d Cumulative Variance',loc),'FontSize',9)
            xlabel('t','FontSize',9)
            hold off
        end

        leg = legend('Standard','BP Knock','Cyp Knock');
        % Programatically move the Legend
        newPosition = [0.9 0.4 0.01 0.1];
        newUnits = 'normalized';
        set(leg,'Position', newPosition,'Units', newUnits,'FontSize',8);

        noteGCF('Cumulative Variance Plots',10,'SPDECumVar', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

    %}


        %% Hox/Krox 2D Gradient Plots
        if enableState3
        figure('Position',[100,100,1480,720])
        subplot(2,3,1)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,HoxGrad)
        colorbar
        title('Gradient Plot, Hox')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        subplot(2,3,2)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,HoxGradBPknock)
        colorbar
        title('Gradient Plot, Hox. BP Knockdown')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        subplot(2,3,3)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,HoxGradCypknock)
        colorbar
        title('Gradient Plot, Hox. Cyp Knockdown')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        subplot(2,3,4)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,KroxGrad)
        colorbar
        title('Gradient Plot, Krox')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        subplot(2,3,5)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,KroxGradBPknock)
        colorbar
        title('Gradient Plot, Krox. BP Knockdown')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        subplot(2,3,6)
        hold on
        imagesc(xStart:dx:xEnd,0:dy:yEnd,KroxGradCypknock)
        colorbar
        title('Gradient Plot, Krox. Cyp Knockdown')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([xStart xEnd 0 yEnd])
        hold off

        noteGCF('End Time Gradients',20,'SPDE2DHKGrad', ...
            id,email,diaryTxt,description,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)
        end
    end
end
%% Mean Var Analysis
if K>0
    for location = 1:4
		margv = .065;
		margh = .045;
		pad = .06;
		map = [0 0 .5; .8 0 0];
		figDefs = get(0,'defaultfigureposition');

		%Zeta Plots
		figure('Position',[figDefs(1),figDefs(2),720,720])
        colormap(map)
        subplot_tight(2,1,1,[margv margh]);
        h1 = histogram(ZetasBP(location,:),10);
        h1.BinWidth = .10;
        hold on
        if knockCyp
            h2 =histogram(ZetasCyp(location,:),10);
            h2.BinWidth = .10;
        end
        
        axis([0 1 0 inf])
        title(sprintf('RA_{in}, X=%d',(location-1)*100),'FontSize',8)
        y = ylabel('Frequency','FontSize',8);
        set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
        x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
        set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
        text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
        legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')

        subplot_tight(2,1,2,[margv margh]);
        h1 = histogram(ZetasBPRAR(location,:),10);
        h1.BinWidth = .10;
        hold on
        if knockCyp
            h2 =histogram(ZetasCypRAR(location,:),10);
            h2.BinWidth = .10;
        end
        axis([0 1 0 inf])
        title(sprintf('RAR, X=%d',(location-1)*100),'FontSize',8)
        y = ylabel('Frequency','FontSize',8);
        set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
        x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
        set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
        text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
        legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')

		pubGCF('',10,sprintf('SPDEKnockZeta%d',location),id,email,diaryTxt,description,location,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

        %{
		%Scatter Plots
		figure('Position',[figDefs(1),figDefs(2),960,720])
		margv = .07;
		margh = .043;
		pad = .06;
		colormap(map)
		hold on
		subplot_tight(2,2,1,[margv margh]);
		box on
		gscatter(DeltaMeansBP(location,:),DeltaVarsBP(location,:),varBPDir(location,:),map)
		axis([0 100 0 100])
		title(sprintf('RA_{in}: BP Knockdown, X=%d',(location-1)*100),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
		legend('Mean Decrease','Mean Increase')
		text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off
        
        if knockCyp
            hold on
            subplot_tight(2,2,2,[margv margh]);
            gscatter(DeltaMeansCyp(location,:),DeltaVarsCyp(location,:),meanCypDir(location,:),map)
            axis([0 100 0 100])
            title(sprintf('RA_{in}: Cyp Knockdown, X=%d',(location-1)*100),'FontSize',8)
            y = ylabel('%\Delta Variance','FontSize',8);
            x = xlabel('%\Delta Mean','FontSize',8);
            set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
            set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
            legend('Mean Decrease','Mean Increase')
            text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
            hold off
        end

		hold on
		subplot_tight(2,2,3,[margv margh]);
		colormap(map)
		gscatter(DeltaMeansBPRAR(location,:),DeltaVarsBPRAR(location,:),varBPDirR(location,:),map)
		axis([0 100 0 100])
		title(sprintf('RAR: BP Knockdown, X=%d',(location-1)*100),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		legend('Var Decrease','Var Increase')
		text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

        if knockCyp
            subplot_tight(2,2,4,[margv margh]);
            gscatter(DeltaMeansCypRAR(location,:),DeltaVarsCypRAR(location,:),meanCypDirR(location,:),map)
            hold on
            axis([0 100 0 100])
            title(sprintf('RAR: Cyp Knockdown, X=%d',(location-1)*100),'FontSize',8)
            y = ylabel('%\Delta Variance','FontSize',8);
            x = xlabel('%\Delta Mean','FontSize',8);
            set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
            set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
            legend('Mean Decrease','Mean Increase')
            text(-.05,1.1,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
            hold off
        end
		pubGCF('',10,sprintf('SPDEKnockScat%d',location),id,email,diaryTxt,description,location,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)
        %}
    end

    %% Boundary Sharpening
    for location=1:length(sharpLocations)
    	margv = .065;
		margh = .045;
		pad = .06;
		map = [0 0 .5; .8 0 0];
        %{
		figDefs = get(0,'defaultfigureposition');
		%Zeta Plots
		figure('Position',[figDefs(1),figDefs(2),960,720])
		hold on
		colormap(map)
		subplot_tight(2,2,1,[margv margh]);
		hist(BZetasBP(location,:))
		axis([0 1 0 inf])
		title(sprintf('RA_{in} Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('Frequency','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')

		subplot_tight(2,2,2,[margv margh]);
		hist(BZetasCyp(location,:))
		axis([0 1 0 inf])
		title(sprintf('RA_{in} Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',4)
		y = ylabel('Frequency','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')

		subplot_tight(2,2,3,[margv margh]);
		hist(BZetasBPRAR(location,:))
		axis([0 1 0 inf])
		title(sprintf('RAR Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('Frequency','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')

		subplot_tight(2,2,4,[margv margh]);
		hist(BZetasCypRAR(location,:))
		axis([0 1 0 inf])
		title(sprintf('RAR Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',4)
		y = ylabel('Frequency','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		x = xlabel('\zeta','FontSize',100,'FontWeight','bold');
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		text(-.05,1.1,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

		pubGCF('',10,sprintf('SPDEBoundaryZetaX%d',location),id,email,diaryTxt,description,location,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

		%Scatter Plots
		figure('Position',[figDefs(1),figDefs(2),960,720])
		margv = .07;
		margh = .043;
		pad = .06;
		colormap(map)
		hold on
		subplot_tight(2,2,1,[margv margh]);
		box on
		gscatter(BDeltaMeansBP(location,:),BDeltaVarsBP(location,:),BvarBPDir(location,:),map)
		axis([0 100 0 100])
		title(sprintf('RA_{in} Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
		legend('Mean Decrease','Mean Increase')
		text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

		hold on
		subplot_tight(2,2,2,[margv margh]);
		gscatter(BDeltaMeansCyp(location,:),BDeltaVarsCyp(location,:),BmeanCypDir(location,:),map)
		axis([0 100 0 100])
		title(sprintf('RA_{in} Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		legend('Mean Decrease','Mean Increase')
		text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

		hold on
		subplot_tight(2,2,3,[margv margh]);
		colormap(map)
		gscatter(BDeltaMeansBPRAR(location,:),BDeltaVarsBPRAR(location,:),BmeanBPDirR(location,:),map)
		axis([0 100 0 100])
		title(sprintf('RAR Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		legend('Var Decrease','Var Increase')
		text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

		subplot_tight(2,2,4,[margv margh]);
		gscatter(BDeltaMeansCypRAR(location,:),BDeltaVarsCypRAR(location,:),BmeanCypDirR(location,:),map)
		hold on
		axis([0 100 0 100])
		title(sprintf('RAR Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('%\Delta Mean','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		legend('Mean Decrease','Mean Increase')
		text(-.05,1.1,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

		pubGCF('',10,sprintf('SPDEBoundaryScatPerX%d',location),id,email,diaryTxt,description,location,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)
        %}
        %Absolte X Scatter Plots
		figure('Position',[figDefs(1),figDefs(2),960,720])
		margv = .07;
		margh = .043;
		pad = .06;
		colormap(map)
		hold on
		subplot_tight(2,2,1,[margv margh]);
		box on
		gscatter(BDeltaMeansBPX(location,:),BDeltaVarsBP(location,:),BvarBPDir(location,:),map)
		ylim([0 100])
        xlim([0 150])
		title(sprintf('RA_{in} Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('\Delta X','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.035, 0]);
		legend('Var Decrease','Var Increase')
		text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

        if knockCyp
            hold on
            subplot_tight(2,2,2,[margv margh]);
            gscatter(BDeltaMeansCypX(location,:),BDeltaVarsCyp(location,:),BvarCypDir(location,:),map)
            ylim([0 100])
            xlim([0 150])
            title(sprintf('RA_{in} Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',8)
            y = ylabel('%\Delta Variance','FontSize',8);
            x = xlabel('\Delta X','FontSize',8);
            set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
            set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
            legend('Var Decrease','Var Increase')
            text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
            hold off
        end
		hold on
		subplot_tight(2,2,3,[margv margh]);
		colormap(map)
		gscatter(BDeltaMeansBPXRAR(location,:),BDeltaVarsBPRAR(location,:),BvarBPDirR(location,:),map)
        ylim([0 100])
        xlim([0 150])
		title(sprintf('RAR Location %.2f: BP Knockdown',sharpLocations(location)),'FontSize',8)
		y = ylabel('%\Delta Variance','FontSize',8);
		x = xlabel('\Delta X','FontSize',8);
		set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
		set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
		legend('Var Decrease','Var Increase')
		text(-.05,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
		hold off

        if knockCyp
            subplot_tight(2,2,4,[margv margh]);
            gscatter(BDeltaMeansCypXRAR(location,:),BDeltaVarsCypRAR(location,:),BvarCypDirR(location,:),map)
            hold on
            ylim([0 100])
            xlim([0 150])
            title(sprintf('RAR Location %.2f: Cyp Knockdown',sharpLocations(location)),'FontSize',8)
            y = ylabel('%\Delta Variance','FontSize',8);
            x = xlabel('\Delta X','FontSize',8);
            set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
            set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
            legend('Var Decrease','Var Increase')
            text(-.05,1.1,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
            hold off
        end

		pubGCF('',10,sprintf('SPDEBoundaryScatX%d',location),id,email,diaryTxt,description,location,tEnd,dt,tTot2,dt2,K,sigma,reaction,regime,noiseRegime,diffOn,vbpAdjust,enableState3)

    end
end

%% Finish

if(batch)
    close all
    delete(gcp)
end

end
