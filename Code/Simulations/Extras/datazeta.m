wt = [0.5783402	0.4911179	0.5546541
0.5428075	0.552307	0.5504302
0.6403488	0.533347	0.476655
0.531559	0.5670277	0.499975
0.571288	0.5309062	0.5919694
0.5245949	0.5842766	0.5742278
0.5722372	0.5581799	0.5334885
0.4742986	0.4958641	0.6099987
0.5205743	0.4567947	0.500144
0.512769	0.5417599	0.4799038
0.5246117	0.6176586	0.5873936
0.4793187	0.562873	0.4984883
0.5060365	0.6199829	0.5836432
0.5555048	0.5872192	0.5545676
0.6075143	0.5480298	0.5347216
0.5623472	0.4979841	0.6465824
0.6187901	0.5497501	0.5522132
0.5337254	0.599875	0.5622866
0.5107802	0.6147456	0.5910501
0.6104913	0.5292761	0.6365815];

crabp_MO = [0.5389089	0.6455694	0.4914
0.4253175	0.5335528	0.5746884
0.7384794	0.4494058	0.5820897
0.5221025	0.419495	0.4441727
0.5237879	0.5743605	0.4311767
0.5614741	0.5165489	0.594525
0.6773779	0.4952258	0.6717257
0.4845487	0.5909755	0.5366718
0.5215672	0.5603854	0.5529989
0.4845109	0.5420983	0.6616215
0.5064151	0.6064063	0.6639866
0.5384201	0.6301973	0.5673419
0.5969467	0.6408162	0.6466376
0.5285755	0.5291795	0.4859329
0.6401639	0.5401555	0.5210929
0.4767209	0.4643835	0.6396599
0.5027425	0.424687	0.5526467
0.5349731	0.4689134	0.4817443
0.5727659	0.5934474	0.5131068
0.6368086	0.487999	0.5277109];

crabp_GOF = [0.5422791	0.5714271	0.5319704
0.5474566	0.5672395	0.5093702
0.5503271	0.535015	0.4740524
0.5383675	0.5054086	0.5894117
0.529705	0.5350677	0.5305242
0.5449305	0.5117889	0.5971853
0.5094531	0.5801013	0.5405144
0.5337871	0.5104733	0.5162523
0.5161992	0.5147779	0.578106
0.5751932	0.5132552	0.4891806
0.5581783	0.556061	0.6430473
0.5635384	0.4689428	0.575885
0.5518768	0.5052181	0.5970805
0.5450082	0.5195428	0.543272
0.6078966	0.5286955	0.5150045
0.5573539	0.4976124	0.5848044
0.5790219	0.5475374	0.5496631
0.5398723	0.4831792	0.4769582
0.5836176	0.5236659	0.5452471
0.5254419	0.6050276	0.5629143];

cyp_MO = [0.4979462	0.6025833	0.5924292
0.6763057	0.6235694	0.6259145
0.5893201	0.5202323	0.5895125
0.6368237	0.7045614	0.5329615
0.6171671	0.5736457	0.53755
0.5681459	0.6989186	0.5357372
0.5271311	0.6033238	0.6149325
0.6782889	0.4861975	0.5720689
0.6224659	0.645	0.6710623
0.5721823	0.5953903	0.6360976
0.6581916	0.4939723	0.6228608
0.5637134	0.6689985	0.6122565
0.6286107	0.5883861	0.6211333
0.6521581	0.5338959	0.5747706
0.6066859	0.5893895	0.5513119
0.6815836	0.6251027	0.5333824
0.6030113	0.640182	0.634354
0.6569985	0.6948902	0.5436
0.6007996	0.548104	0.6345961
0.5363169	0.5733831	0.7171537];

cyp_GOF = [0.4853933	0.4208688	0.5408911
0.5028391	0.5228376	0.5224542
0.444964	0.4363813	0.492255
0.4303466	0.4437294	0.4312733
0.5381461	0.4826587	0.5307346
0.5424391	0.5332131	0.5094891
0.4776112	0.4804906	0.4076749
0.5097582	0.5477998	0.4419229
0.4867272	0.5157996	0.4786709
0.4948456	0.5304875	0.5134684
0.5482921	0.444629	0.4796876
0.4935842	0.439197	0.4665237
0.4239097	0.5644578	0.5505626
0.5233691	0.4747096	0.554074
0.5138848	0.4720975	0.5846158
0.4334305	0.442798	0.5292293
0.4060438	0.4073365	0.4513608
0.5511927	0.4960784	0.4703161
0.5973376	0.5278941	0.492452
0.5725319	0.4453091	0.5307257];

zeroOffset = 0;
wt = wt - zeroOffset;
crabp_MO = crabp_MO - zeroOffset;
crabp_GOF  = crabp_GOF - zeroOffset;
cyp_MO = cyp_MO - zeroOffset;
cyp_GOF = cyp_GOF -zeroOffset;

% Pooled analysis

fprintf('Pooled analysis\n')

mean1 = mean(wt(:));
var1  = var(wt(:));
cov1 = sqrt(var1)/mean1;

mean21 = mean(crabp_MO(:));
var21  = var(crabp_MO(:));
mean22 = mean(crabp_GOF(:));
var22  = var(crabp_GOF(:));

mean31 = mean(cyp_MO(:));
var31  = var(cyp_MO(:));
mean32 = mean(cyp_GOF(:));
var32  = var(cyp_GOF(:));
cov31 = sqrt(var31)/mean31;
cov32 = sqrt(var32)/mean32;

perChangeMeanBP21 = abs((mean21-mean1)./max(mean1,mean21) * 100);
perChangeVarBP21  = abs((var21-var1)./max(var1,var21) * 100);
ZetasBP1 =perChangeVarBP21./(perChangeMeanBP21+perChangeVarBP21);
perChangeMeanBP22 = abs((mean22-mean1)./max(mean1,mean22) * 100);
perChangeVarBP22  = abs((var22-var1)./max(var1,var22) * 100);
ZetasBP2 =perChangeVarBP22./(perChangeMeanBP22+perChangeVarBP22);
ZetaBP = [ZetasBP1 ZetasBP2]

perChangeMeanCyp31 = abs((mean31-mean1)./max(mean1,mean31) * 100);
perChangeVarCyp31  = abs((var31-var1)./max(var1,var31) * 100);
perChangeCovCyp31  = abs((cov31-cov1)./max(cov1,cov31) * 100);
ZetasCyp1 =perChangeVarCyp31./(perChangeMeanCyp31+perChangeVarCyp31);
ZetasCovCyp1 =perChangeCovCyp31./(perChangeMeanCyp31+perChangeCovCyp31);
perChangeMeanCyp32 = abs((mean32-mean1)./max(mean1,mean32) * 100);
perChangeVarCyp32  = abs((var32-var1)./max(var1,var32) * 100);
perChangeCovCyp32  = abs((cov32-cov1)./max(cov1,cov32) * 100);
ZetasCyp2 =perChangeVarCyp32./(perChangeMeanCyp32+perChangeVarCyp32);
ZetasCovCyp2 =perChangeCovCyp32./(perChangeMeanCyp32+perChangeCovCyp32);
ZetaCyp = [ZetasCyp1 ZetasCyp2]
ZetaCovCyp =  [ZetasCovCyp1 ZetasCovCyp2]

% Individual

fprintf('Individual analysis\n')

mean1 = mean(wt);
var1  = var(wt);
cov1 = sqrt(var1)./mean1;

mean1 = [mean1(1) mean1(1) mean1(1) mean1(2) mean1(2) mean1(2) mean1(3) mean1(3) mean1(3)];
var1 = [var1(1) var1(1) var1(1) var1(2) var1(2) var1(2) var1(3) var1(3) var1(3)];

mean21 = mean(crabp_MO);
var21  = var(crabp_MO);
mean22 = mean(crabp_GOF);
var22  = var(crabp_GOF);

mean21 = [mean21(1) mean21(2) mean21(3) mean21(1) mean21(2) mean21(3) mean21(1) mean21(2) mean21(3)];
mean22 = [mean22(1) mean22(2) mean22(3) mean22(1) mean22(2) mean22(3) mean22(1) mean22(2) mean22(3)];
mean222= [mean22(1) mean22(1) mean22(1) mean22(2) mean22(2) mean22(2) mean22(3) mean22(3) mean22(3)];

var21 = [var21(1) var21(2) var21(3) var21(1) var21(2) var21(3) var21(1) var21(2) var21(3)];
var22 = [var22(1) var22(2) var22(3) var22(1) var22(2) var22(3) var22(1) var22(2) var22(3)];
var222 = [var22(1) var22(1) var22(1) var22(2) var22(2) var22(2) var22(3) var22(3) var22(3)];

mean31 = mean(cyp_MO);
var31  = var(cyp_MO);
mean32 = mean(cyp_GOF);
var32  = var(cyp_GOF);
cov31 = sqrt(var31)./mean31;
cov32 = sqrt(var32)./mean32;

mean31 = [mean31(1) mean31(2) mean31(3) mean31(1) mean31(2) mean31(3) mean31(1) mean31(2) mean31(3)];
mean32 = [mean32(1) mean32(2) mean32(3) mean32(1) mean32(2) mean32(3) mean32(1) mean32(2) mean32(3)];
mean322 = [mean32(1) mean32(1) mean32(1) mean32(2) mean32(2) mean32(2) mean32(3) mean32(3) mean32(3)];

var31 = [var31(1) var31(2) var31(3) var31(1) var31(2) var31(3) var31(1) var31(2) var31(3)];
var32 = [var32(1) var32(2) var32(3) var32(1) var32(2) var32(3) var32(1) var32(2) var32(3)];
var322 = [var32(1) var32(1) var32(1) var32(2) var32(2) var32(2) var32(3) var32(3) var32(3)];

perChangeMeanBP21 = abs((mean21-mean1)./max(mean1,mean21) * 100);
perChangeVarBP21  = abs((var21-var1)./max(var1,var21) * 100);
ZetasBP1 =perChangeVarBP21./(perChangeMeanBP21+perChangeVarBP21);
perChangeMeanBP22 = abs((mean22-mean1)./max(mean1,mean22) * 100);
perChangeVarBP22  = abs((var22-var1)./max(var1,var22) * 100);
ZetasBP2 =perChangeVarBP22./(perChangeMeanBP22+perChangeVarBP22);
perChangeMeanBP222 = abs((mean222-mean21)./max(mean21,mean22) * 100);
perChangeVarBP222  = abs((var222-var21)./max(var1,var22) * 100);
ZetasBP22 =perChangeVarBP222./(perChangeMeanBP222+perChangeVarBP222);
ZetaBP = [ZetasBP1 ZetasBP2 ZetasBP22]

fprintf('Best ZetaBP estimate:%f\n',mean(ZetaBP))
perChangeMeanCyp31 = abs((mean31-mean1)./max(mean1,mean31) * 100);
perChangeMeanCyp312 = abs((mean31-mean322)./max(mean322,mean31) * 100);
perChangeVarCyp31  = abs((var31-var1)./max(var1,var31) * 100);
perChangeVarCyp312  = abs((var31-var322)./max(var322,var31) * 100);
%perChangeCovCyp31  = abs((cov31-cov1)./max(cov1,cov31) * 100);
ZetasCyp1 =perChangeVarCyp31./(perChangeMeanCyp31+perChangeVarCyp31);
ZetasCyp12 =perChangeVarCyp312./(perChangeMeanCyp312+perChangeVarCyp312);
%ZetasCovCyp1 =perChangeCovCyp31./(perChangeMeanCyp31+perChangeCovCyp31);
perChangeMeanCyp32 = abs((mean32-mean1)./max(mean1,mean32) * 100);
perChangeVarCyp32  = abs((var32-var1)./max(var1,var32) * 100);
%perChangeCovCyp32  = abs((cov32-cov1)./max(cov1,cov32) * 100);
ZetasCyp2 =perChangeVarCyp32./(perChangeMeanCyp32+perChangeVarCyp32);
%ZetasCovCyp2 =perChangeCovCyp32./(perChangeMeanCyp32+perChangeCovCyp32);
ZetaCyp = [ZetasCyp1 ZetasCyp2 ZetasCyp12]
fprintf('Best ZetaCyp estimate:%f\n',mean(ZetaCyp))
%ZetaCovCyp =  [ZetasCovCyp1 ZetasCovCyp2]

%% Scatter Plots

figure('Position',[figDefs(1),figDefs(2),1000,720])
subplot(2,2,1)
scatter(mean1(1:3:end),var1(1:3:end),100,'filled')
hold on
scatter(mean31(1:3),var31(1:3),100,'filled')
scatter(mean32(1:3),var32(1:3),100,'filled')
y = ylabel(sprintf('Variance Relative RA Abundance',fntSze));
set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'FontWeight','bold');
x = xlabel(sprintf('Mean of Relative RA Abundance',zetaSze),'FontWeight','bold');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
legend({'Wildtype','Knockdown','GoF'},'location','northwest')
hold off


subplot(2,2,2)
scatter(perChangeMeanCyp31,perChangeVarCyp31,100,'filled')
hold on
axis([0 100 0 100])
scatter(perChangeMeanCyp32,perChangeVarCyp32,100,'filled')
scatter(perChangeMeanCyp312,perChangeVarCyp312,100,'filled')
y = ylabel(sprintf('%% Change Variance Relative RA Abundance',fntSze));
set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'FontWeight','bold');
x = xlabel(sprintf('%% Change Mean of Relative RA Abundance',zetaSze),'FontWeight','bold');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
legend({'Wildtype vs Knockdown','Wildtype vs GoF','Knockdown vs GoF'})
hold off

%% STD

std1 = sqrt(var1);
std31 = sqrt(var31);
std32 = sqrt(var32);
std322 = sqrt(var322);
perChangeStdCyp31  = abs((std31-std1)./max(std1,std31) * 100);
perChangeStdCyp312  = abs((std31-std322)./max(std322,std31) * 100);
perChangeStdCyp32  = abs((std32-std1)./max(std1,std32) * 100);

figure('Position',[figDefs(1),figDefs(2),1000,720])
subplot(2,2,1)
scatter(mean1(1:3:end),std1(1:3:end),100,'filled')
hold on
scatter(mean31(1:3),std31(1:3),100,'filled')
scatter(mean32(1:3),std32(1:3),100,'filled')
y = ylabel(sprintf('Standard Deviation Relative RA Abundance',fntSze));
set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'FontWeight','bold');
x = xlabel(sprintf('Mean of Relative RA Abundance',zetaSze),'FontWeight','bold');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
legend({'Wildtype','Knockdown','GoF'},'location','northwest')
hold off


subplot(2,2,2)
scatter(perChangeMeanCyp31,perChangeStdCyp31,100,'filled')
hold on
axis([0 100 0 100])
scatter(perChangeMeanCyp32,perChangeStdCyp312,100,'filled')
scatter(perChangeMeanCyp312,perChangeStdCyp32,100,'filled')
y = ylabel(sprintf('%% Change Standard Deviation Relative RA Abundance',fntSze));
set(y, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'FontWeight','bold');
x = xlabel(sprintf('%% Change Mean of Relative RA Abundance',zetaSze),'FontWeight','bold');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
legend({'Wildtype vs Knockdown','Wildtype vs GoF','Knockdown vs GoF'})
hold off

%% Load Data

addRA     = load('D:\OneDrive\Important_Projects\MeanVar Independance\arxData\Hists\16-08-13-225901-5-RAsde.mat');
multRA    = load('D:\OneDrive\Important_Projects\MeanVar Independance\arxData\MultiHists\16-08-15-094454-51-RAsde.mat');
multRAR   = load('D:\OneDrive\Important_Projects\MeanVar Independance\arxData\noiseLocality\16-08-14-124244-52-RAsde.mat');
multRAout = load('D:\OneDrive\Important_Projects\MeanVar Independance\arxData\noiseLocality\16-10-26-125147-58-RAsde.mat');
multRABP = load('D:\OneDrive\Important_Projects\MeanVar Independance\arxData\noiseLocality\16-10-26-134156-54-RAsde.mat');

data = {addRA,multRA,multRAR,multRAout,multRABP};
descriptions = {'Additive on RA','Multiplicative on RA','Multiplicative on RAR','Multiplicative on RAout','Multiplicative on RABP','Experimental Values'};
%% Plot Parameters

id = timeStamp();
figDefs = get(0,'defaultfigureposition');
margvHist = .08;
marghHist = .02;
margvHistBot = .065;
padHist = .06;
margvScat = .07;
marghScat = .08;
padScat = .06;
fntSze = 18;
titSze = 15;
lw = 3;
zetaSze = 25;
letterSize = 12;
histBins = 10;
map = [0 0 .5; .8 0 0];
varginTxt = false;
extendedScatter = false;
edges = 0:.1:1;
N = 1000;
bw = .1;
ZetaCyp = [0.7049    0.8506    0.6770    0.7461    0.8644    0.7364    0.6244    0.8525    0.5112    0.7090    0.0325    0.0958    0.7316    0.3024 0.3053    0.5591    0.4324    0.5273    0.0808    0.6696    0.2718    0.5254    0.7107    0.4630    0.5730    0.7449    0.5201];

%%

subplot(2,2,3)
hold on
maxf = 0;
for i = 1:length(data)
    [f,xi] = ksdensity(data{i}.ZetasCyp,'support','positive','width',bw);
    maxf = max(max(maxf,f));
    semilogy(xi,f,'LineWidth',3)
    y = ylabel(sprintf('Estimated Probability Distribution Function for \\fontsize{%d}\\zeta',fntSze));
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',fntSze));
    axis([0 1 0 5])
    set(gca,'FontSize',fntSze)
end
y1=get(gca,'ylim');
%[f,xi] = ksdensity(ZetaCyp,'support','positive','width',0.1);
%semilogy(xi,f,'r--','LineWidth',3)
for i = 1:length(ZetaCyp)
    SP=ZetaCyp(i); %your point goes here 
    plot([0 SP],[0 0],'ro','LineWidth',3);
end
legend(descriptions)
hold off

subplot(2,2,4)
hold on
maxf = 0;
for i = 1:length(data)
    [f,xi] = ksdensity(data{i}.ZetasCyp,'support','positive','width',bw,'function','cdf');
    maxf = max(max(maxf,f));
    semilogy(xi,f,'LineWidth',3)
    y = ylabel(sprintf('Estimated Cumulative Distribution for \\fontsize{%d}\\zeta',fntSze));
    x = xlabel(sprintf('\\fontsize{%d}\\zeta',fntSze));
    axis([0 1 0 1])
    set(gca,'FontSize',fntSze)
end
[f,xi] = ksdensity(ZetaCyp,'support','positive','width',0.1,'function','cdf');
maxf = max(max(maxf,f));
semilogy(xi,f,'r--','LineWidth',3)
legend(descriptions)
y1=get(gca,'ylim');
hold off

%{
figure
maxf = 0;
for i = 1:length(data)
    hold on
    [f,xi] = ksdensity(data{i}.ZetasCyp2,'support','positive','width',bw);
    maxf = max(max(maxf,f));
    semilogx(xi,f,'LineWidth',3)
    y = ylabel(sprintf('Probability of \\fontsize{%d}\\zeta_{CoV}',fntSze));
    x = xlabel(sprintf('\\fontsize{%d}\\zeta_{CoV}',fntSze));
    axis([0 1 0 inf])
    set(gca,'FontSize',fntSze)
end

legend(descriptions)
SP=.293; %your point goes here 
plot([SP SP],[0 maxf],'r--','LineWidth',3);
hold off
%}
