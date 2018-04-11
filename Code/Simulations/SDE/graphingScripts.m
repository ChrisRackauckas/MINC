    margv = .065;
    margh = .045;
    pad = .06;
    
    %Zeta Plots
    map = [0 0 .5; .8 0 0];
    
    figDefs = get(0,'defaultfigureposition');
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(map)
    subplot_tight(2,1,1,[margv margh]);
    h1 = histogram(ZetasBP,10)
    hold on 
    h2 =histogram(ZetasCyp,10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Initial State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')
    
    figDefs = get(0,'defaultfigureposition');
    subplot_tight(2,1,2,[margv margh]);
    h1 = histogram(ZetasBP,10)
    hold on 
    h2 =histogram(ZetasCyp,10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Terminal State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')
    
    
    
    
    
    
    
    
    
    
    
    figDefs = get(0,'defaultfigureposition');
    figure('Position',[figDefs(1),figDefs(2),720*2,720])
    hold on
    subaxis(1,2,1, 'Spacing', 0.00, 'Padding', .05, 'Margin', 0.01);
    colormap(map)
    gscatter(DeltaMeansCyp,DeltaVarsCyp,meanCypDir,map)
    axis([0 100 0 100])
    title('RA_{in}: Cyp Knockdown','FontSize',8)
    ylabel('%\Delta Variance','FontSize',8)
    xlabel('%\Delta Mean','FontSize',8)
    legend('Mean Decrease','Mean Increase')
    text(-.115,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    
    hold on
    subaxis(1,2,2, 'Spacing', 0.00, 'Padding', .05, 'Margin', 0.01);
    colormap(map)
    gscatter(DeltaMeansCyp,DeltaVarsCyp,meanCypDir,map)
    axis([0 100 0 100])
    title('RA_{in}: Cyp Knockdown','FontSize',8)
    ylabel('%\Delta Variance','FontSize',8)
    xlabel('%\Delta Mean','FontSize',8)
    legend('Mean Decrease','Mean Increase')
    text(-.115,1.1,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    
    
    %% fig 4 5
    
        margv = .065;
    margh = .045;
    pad = .06;
    
    %Zeta Plots
    map = [0 0 .5; .8 0 0];
    
    figDefs = get(0,'defaultfigureposition');
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(map)
    subplot_tight(2,1,1,[margv margh]);
    h1 = histogram(ZetasBPRAR,10)
    hold on 
    h2 =histogram(ZetasCypRAR,10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Initial State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')
    
    figDefs = get(0,'defaultfigureposition');
    subplot_tight(2,1,2,[margv margh]);
    h1 = histogram(ZetasBPRAR,10)
    hold on 
    h2 =histogram(ZetasCypRAR,10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Terminal State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')
    
    
    
    %% SPDE
    
    
        margv = .065;
    margh = .045;
    pad = .06;
    
    %Zeta Plots
    map = [0 0 .5; .8 0 0];
    
    figDefs = get(0,'defaultfigureposition');
    figure('Position',[figDefs(1),figDefs(2),720,720])
    
    colormap(map)
    subplot_tight(2,1,1,[margv margh]);
    h1 = histogram(ZetasBP(1,:),10)
    hold on 
    h2 =histogram(ZetasCyp(1,:),10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Initial State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')
    
    figDefs = get(0,'defaultfigureposition');
    subplot_tight(2,1,2,[margv margh]);
    h1 = histogram(ZetasBP(4,:),10)
    hold on 
    h2 =histogram(ZetasCyp(4,:),10)
    h1.BinWidth = .10;
    h2.BinWidth = .10;
    axis([0 1 0 inf])
    title('Terminal State: BP Knockdown','FontSize',8)
    y = ylabel('Frequency','FontSize',8);
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    x = xlabel('\zeta','FontSize',8,'FontWeight','bold');
    set(x, 'Units', 'Normalized', 'Position', [0.5, -0.025, 0]);
    text(-.05,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12,'fontw','b')
    
    legend('Cascade-Strength Knockdown','Feedback-Strength Knockdown')