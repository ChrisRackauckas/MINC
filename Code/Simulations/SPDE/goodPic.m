   varmult = 125;
    vars1 = [   0.000849327635420
   0.168467447190003
   varmult*0.100000000000000
  49.788238866358910
   0.012505206316120
   0.001403947874187
   0.000004323632946
   0.759871222501553
   0.000615282202737
   4.473519914922131
   0.103662576283866
   0.084899283524205
   0.007407926966825
   0.000251061485367
   0.116941937282287
   0.000029174670369
   0.019843436847897
   0.009805651879419
   varmult*0.001000000000000
   varmult*0.001000000000000];


    vars2 = [   0.001348528850545
   0.002196437694883
   varmult*0.100000000000000
   0.004070277191119
   0.036845147022038
   0.000123587462458
   0.000195235808878
   0.001384785085646
   0.091600171688902
   2.642862950864342
   0.000629326548045
   0.033529481558713
   0.000281011636613
   0.000031310938709
   0.000005670662321
   0.000681070076220
   0.000007821564772
   0.000049249468243
   varmult*0.001000000000000
   varmult*0.001000000000000];
    
    vars3 = [  0.004525750147322
   0.004272950914649
   varmult*0.100000000000000
   9.529657183225945
   0.019221199393119
   0.000519424103233
   0.016085607696992
   1.586466742119592
   0.000023208008693
   0.017968375922538
   0.001560517651868
   0.000005671785089
   0.009585351430661
   0.025194218801073
   0.000248953561258
   0.000242109422453
   0.000008924817407
   0.211208146700360
   varmult*0.001000000000000
   varmult*0.001000000000000];

    vars4 = [   0.066067259198318
   0.001369222751244
   varmult*0.100000000000000
   0.015776092522099
   0.003562529255699
   0.001702100058829
   0.031921820900823
   0.116619377320080
   0.000076341211651
   0.410227667148109
   0.002053293373190
   0.029400997263032
   0.168131815498571
   0.000132388966949
   0.000033331716622
   0.108522584432902
   0.006354576100516
   0.035375418104496
   varmult*0.001000000000000
   varmult*0.001000000000000];

    vars5 = [   0.295129189652166
   0.000464268965837
   varmult*0.100000000000000
   3.860612935399694
   0.020582708360037
   0.000055080633863
   0.119020415179069
   0.013874782200846
   0.000535099405926
   0.516475547476636
   0.213299925415754
   0.000267562388833
   0.028482301823693
   0.005883447633690
   0.011735500099824
   0.004637374377356
   0.259592087044417
   0.000012160701113
   varmult*0.001000000000000
   varmult*0.001000000000000];
    
    vars6 =1.0e+02 * [0.000202426760945
   0.001252063004959
   varmult*0.001000000000000
   2.786101405164940
   0.000102310592485
   0.000000434537191
   0.000000065279001
   0.000826204017668
   0.000000366191550
   0.033878800806565
   0.000508605594200
   0.000005198898552
   0.002171424767104
   0.000001731681665
   0.000002056057290
   0.003132506307883
   0.000055942972016
   0.000009340849798
   varmult*0.000010000000000
   varmult*0.000010000000000];

%%

[RAout,RA,RABP,RAR,R,BP,inits] = simMotif52(vars1,5e-5,2)

%%


p = convertVarsToP(vars1);

[~,~,RAinGrad,RARGrad,RAinGradBPknock,RARGradBPknock] = driver('Test Run',2e6,5e-5,0,0,'MultIn',0,0,'diffOn',true,'verbose',true,'initRun',true,'knock',0.0001,'suppliedP',p);

        %% 2D Gradient Plots
        fntSze = 14;
        XimgStart = 200;
        XimgEnd = 330;
        xStart = -100;
        xEnd = 400;
        dx = 5;
        dy = 5;
        yEnd = 100;

        StartIdx = (XimgStart-xStart)/dx;
        EndIdx = length(xStart:dx:xEnd)-(xEnd-XimgEnd)/dx;
        
        figure('Position',[100,100,1480,720])
        subplot(2,2,1)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAinGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAin')
        xlabel('x','Fontsize',1)
        ylabel('y','Fontsize',1)
        axis([XimgStart XimgEnd 0 yEnd])
        set(gca,'FontSize',fntSze)
        %caxis([.052  .065])
        hold off

        subplot(2,2,2)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RARGrad(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAR')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        set(gca,'FontSize',fntSze)
        caxis([2e-5 3e-5])
        hold off

        subplot(2,2,3)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RAinGradBPknock(:,StartIdx:EndIdx,:))
        colorbar
        title('Gradient Plot, RAin')
        xlabel('x','Fontsize',1)
        ylabel('y','Fontsize',1)
        axis([XimgStart XimgEnd 0 yEnd])
        set(gca,'FontSize',fntSze)
        caxis([.052 .065])
        hold off

        subplot(2,2,4)
        hold on
        imagesc(XimgStart:dx:XimgEnd,0:dy:yEnd,RARGradBPknock(:,StartIdx:EndIdx))
        colorbar
        title('Gradient Plot, RAR')
        xlabel('x','Fontsize',8)
        ylabel('y','Fontsize',8)
        axis([XimgStart XimgEnd 0 yEnd])
        set(gca,'FontSize',fntSze)
        caxis([2e-5 3e-5])
        hold off

%%

p = struct();
p.prodLimit = 200;
p.epsinMult = 12.5000;
p.epsinAdd = 0;
p.epsoutMult = 2*6.25;
p.epsoutAdd = 0;
p.epsRAdd = 0;
p.epsRMult = 2*6.25;
p.epsHAdd = 0;
p.epsHMult = 0;
p.epsKAdd = 0;
p.epsKMult = 0;
p.DRA = 20;
p.Vbp = 0.3046*100;
p.gamma = 5.7150e-07;
p.beta = 0.0542;
p.kp = 1.3303e-06;
p.kmax = 4.4936e-06;
p.beta0 = 8.4704e-06;
p.lambda = 2.1760e-07;
p.kdeg = 2.1317e-04;
p.vmax = 1.0459e-04;
p.ron = 4.7178e-06;
p.roff = 2.4384e-07;
p.mon = 3.3928e-06;
p.moff = 2.1878e-06;
p.jalpha = 1.7460e-06;
p.jbeta = 1.6315e-07;
p.bpdeg1 = 0.0172;
p.bpdeg2 = 6.9249e-07;
p.Vr = 0.1331;
p.rdeg1 = 0.0073;
p.rdeg2 = 0.0935;
p.d = 8.9290e-07;
p.e = 7.2295e-06;
p.cH = 2.7531e-07;
p.cK = 3.0530e-06;
p.kappaH = 9.9870e-07;
p.kappaK = 1.0394e-07;
p.dH = 3.7919e-06;
p.dK = 1.4014e-07;
p.aH = 3.2693e-07;
p.aK = 9.1342e-06;
p.f0 = 500;

driver('Test Run',2e6,1e-6,0,0,'MultIn',0,0,'diffOn',true,'verbose',true,'initRun',1,'knock',0.01,'suppliedP',p);


%%
driver('Test Run',2e7,1e-6,0,0,'MultIn',0,0,'diffOn',true,'verbose',true,'initRun',true);
