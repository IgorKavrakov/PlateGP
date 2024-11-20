function Example2_Plots(GP,AP)
%% This is a specific function for plots in the paper
close all;

%% Learned flexural rigidity
fprintf('Flexural Rigidity:\n')
fprintf('Learn No BC: (Standard) D/D_true= %2.3f \t\t (MCMC) D/D_true= %2.3f\n',GP(1).Par.D/AP.D, GP(1).Par.D_MCMC/AP.D)
fprintf('Learn Yes BC: (Standard) D/D_true= %2.3f \t\t (MCMC) D/D_true= %2.3f\n',GP(2).Par.D/AP.D, GP(2).Par.D_MCMC/AP.D)

%%
for i = 1:2
    %% Switch between standard and Bayesian optim
    switch i
        case 1 % Standard
            flagMCMC = false;
            predMeshMean = 'f_mesh';
            predMeshStd = 'f_mesh_std';
            subName = '';
        case 2 % MCMC
            flagMCMC = true;
            predMeshMean = 'f_mesh_MCMC';
            predMeshStd = 'f_mesh_std_MCMC';
            subName = '_MCMC';
    end
    %% Analytical solution for ploting countour plots
    AP.x=linspace(0,AP.a,51);AP.y=linspace(0,AP.b,51); %Learning coord
    [AP] = Example2_Analytical_FixedRecPlate(AP);

    %% Figure 1 - Contour Displacements w
    cont_plot(AP.x_mesh,AP.y_mesh,AP.w,...
              GP(1).Train.x_mesh.w,GP(1).Train.y_mesh.w,...
              GP(1).Pred.x_mesh.w,GP(1).Pred.y_mesh.w,GP(1).Pred.(predMeshMean).w,GP(2).Pred.(predMeshMean).w,...
              GP(2).Train.X.w(GP(2).Kernel.BC.w),GP(2).Train.Y.w(GP(2).Kernel.BC.w),...
              [-0.2 1.2],[-0.2:0.2:1.2],...
              '$w/\max |w_{\mathrm{ana}}|$',sprintf('Ex2_Contour_w%s', subName))
          
    %% Figure 2 - Contour Rotation Rx
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Rx,...
              GP(1).Train.x_mesh.Rx,GP(1).Train.y_mesh.Rx,...
              GP(1).Pred.x_mesh.Rx,GP(1).Pred.y_mesh.Rx,GP(1).Pred.(predMeshMean).Rx,GP(2).Pred.(predMeshMean).Rx,...
              GP(2).Train.X.Rx(GP(2).Kernel.BC.Rx),GP(2).Train.Y.Rx(GP(2).Kernel.BC.Rx),...
              [-1.2 1.2],[-1.2:0.4:1.2],...
              '$r_x/\max |r_{x,\mathrm{ana}}|$',sprintf('Ex2_Contour_Rx%s', subName))

    %% Figure 3 - Contour Curvature Kx
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Kx,...
              GP(1).Train.x_mesh.Kx,GP(1).Train.y_mesh.Kx,...
              GP(1).Pred.x_mesh.Kx,GP(1).Pred.y_mesh.Kx,GP(1).Pred.(predMeshMean).Kx,GP(2).Pred.(predMeshMean).Kx,...
              [],[],...
              [-1.2 0.6],[-1.2:0.2:0.6],...
              '$\kappa_x/\max|\kappa_{x,\mathrm{ana}}|$',sprintf('Ex2_Contour_Kx%s', subName))

    %% Figure 4 - Line plot Displacement w 
    line_plot(AP.x,AP.y,AP.w,AP.w,...
              GP(1).Train.x.w,GP(1).Train.y.w,GP(1).Train.f_mesh.w,[],...
              GP(1).Pred.x.w,GP(1).Pred.y.w,GP(1).Pred.(predMeshMean).w,GP(1).Pred.(predMeshStd).w,GP(2).Pred.(predMeshMean).w,GP(2).Pred.(predMeshStd).w,...
              GP(1).Pred.(predMeshMean).w,GP(1).Pred.(predMeshStd).w,GP(2).Pred.(predMeshMean).w,GP(2).Pred.(predMeshStd).w,...
              [0 1],...
              {[0 1;0 0],[GP(2).Train.x.w;GP(2).Train.x.w*0]},...
              '$x/a$ [-]','$w/\max |w_{\mathrm{ana}}|$ [-]','$x/a$ [-]','$w/\max |w_{\mathrm{ana}}|$ [-]',...
              {[0 1],[-0.20 1.20],[0 1],[-0.20 0.20]},{-0.2:0.2:1.2,-0.2:0.1:0.2},...
              {'True','Boundary Conditions (BC)','Observation','Prediction (No BC)','Prediction (Yes BC)'},...
              {'True','Boundary Conditions (BC)','Prediction (No BC)','Prediction (Yes BC)'},...
              sprintf('Ex2_Line_w%s', subName));  
          
    %% Figure 5 - Line plot Rotation Rx 
    line_plot(AP.x,AP.y,AP.Rx,AP.Rx,...
              [],[],[],[],...
              GP(1).Pred.x.Rx,GP(1).Pred.y.Rx,GP(1).Pred.(predMeshMean).Rx,GP(1).Pred.(predMeshStd).Rx,GP(2).Pred.(predMeshMean).Rx,GP(2).Pred.(predMeshStd).Rx,...
              GP(1).Pred.(predMeshMean).Rx,GP(1).Pred.(predMeshStd).Rx,GP(2).Pred.(predMeshMean).Rx,GP(2).Pred.(predMeshStd).Rx,...
              [0 1],...
              {[0 1;0 0],[GP(2).Train.x.w;GP(2).Train.x.w*0]},...
              '$x/a$ [-]','$r_x/\max |r_{x,\mathrm{ana}}|$ [-]','$x/a$ [-]','$r_x/\max |r_{x,\mathrm{ana}}|$ [-]',...
              {[0 1],[-1.20 1.20],[0 1],[-0.20 0.20]},{-1.2:0.4:1.2,-0.2:0.1:0.2},...
              {'True','Boundary Conditions (BC)','Prediction (No BC)','Prediction (Yes BC)'},...
              {'True','Boundary Conditions (BC)','Prediction (No BC)','Prediction (Yes BC)'},...
              sprintf('Ex2_Line_Rx%s', subName));            
          
    %% Figure 6 - Line plot Curvature Kx 
    line_plot(AP.x,AP.y,AP.Kx,AP.Kx,...
              GP(1).Train.x.Kx,GP(1).Train.y.Kx,GP(1).Train.f_mesh.Kx,[],...
              GP(1).Pred.x.Kx,GP(1).Pred.y.Kx,GP(1).Pred.(predMeshMean).Kx,GP(1).Pred.(predMeshStd).Kx,GP(2).Pred.(predMeshMean).Kx,GP(2).Pred.(predMeshStd).Kx,...
              GP(1).Pred.(predMeshMean).Kx,GP(1).Pred.(predMeshStd).Kx,GP(2).Pred.(predMeshMean).Kx,GP(2).Pred.(predMeshStd).Kx,...
              [0 1],...
              {[],[]},...
              '$x/a$ [-]','$\kappa_{x}/\max |\kappa_{x,\mathrm{ana}}|$ [-]','$x/a$ [-]','$\kappa_{x}/\max |\kappa_{x,\mathrm{ana}}|$ [-]',...
              {[0 1],[-1.20 0.61],[0 1],[-0.20 0.20]},{-1.2:0.2:0.6,-0.2:0.1:0.2},...
              {'True','Observation','Prediction (No BC)','Prediction (Yes BC)'},...         
              {'True','Prediction (No BC)','Prediction (Yes BC)'},...         
              sprintf('Ex2_Line_Kx%s', subName));          

    %% Figure 8 - Comparison with & without BCs - only for MCMC
    if flagMCMC
        histo_plot(AP.D,GP(1).Par.D,exp(GP(1).Par.hyp_MCMC(:,4)),...
                        GP(2).Par.D,exp(GP(2).Par.hyp_MCMC(:,4)),...
                sprintf('Ex2_Histogram_D%s',subName))
    end
        
    %% Figure 9 - Correlation matrix, double plots
    if flagMCMC
        correlMat_doubleplot(GP(1).Par.hyp_MCMC(2000:end, 1:4), GP(2).Par.hyp_MCMC(2000:end, 1:4), ...
                {'$A^2$','$l_x$','$l_y$','$D$'}, ...
                sprintf('Ex2_CorrelMat_Double%s', subName))
    end

end

end

%%
function correlMat_doubleplot(dataA, dataB, labels, plotName)
    Colors=[57  106 177
        218 124 48
        62  150 81
        204 37  41
        83  81  84
        107 76  154
        146 36  40
        148 139 61]./255;
    
    numplots = length(labels);
    FigWidth=9; nby=numplots; nbx=numplots; spacey=0.15; spacex=0.15; leftmargin=0.7; rightmargin=0.7; topmargin=0.6; bottommargin=0.6;
    FigDepth=9;
    FontAxis=7; FontLabel=9; FontLegend=7;
    f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
    [ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);
    positions = flip(positions);
    n = 1;
    coefMatA = corrcoef(dataA);
    coefMatB = corrcoef(dataB);
    for i = 1:length(positions)
        for j = 1:length(positions)
            % Current axes
            ax=axes('position',positions{n},'Layer','top'); hold on; box on; grid on;
            % Ticks and labels
            if i == 1 
                xlabel(labels{j}, 'interpreter', 'latex'); 
            end
            if j == 4 
                ylabel(labels{i}, 'interpreter', 'latex'); 
            end
            xticks([])
            yticks([])
            if i == 4; set(gca,'xaxisLocation','top'); xlabel(labels{j}, 'interpreter', 'latex'); end
            if j == 1; set(gca,'yaxisLocation','right'); ylabel(labels{i}, 'interpreter', 'latex'); end
            % Plots
            if i == j % Main diagonal: histogram
                nBins = 50;
                pA = fitdist(dataA(:,i),'Kernel');
                pB = fitdist(dataB(:,i),'Kernel');
                histogram(dataA(:,i), 'numbins', nBins, 'edgecolor', 'none', 'facecolor' ,Colors(1,:), 'facealpha', 0.25, 'Normalization','pdf');
                histogram(dataB(:,i), 'numbins', nBins, 'edgecolor', 'none', 'facecolor' ,Colors(4,:), 'facealpha', 0.25, 'Normalization','pdf');
                x = linspace(ax.XLim(1), ax.XLim(2), 1e3);
                yLim = ax.YLim;
                yA = pdf(pA,x);
                yB = pdf(pB,x);
                plot(x, yA, 'Color', Colors(1,:));
                plot(x, yB, 'Color', Colors(4,:));
                ylim(yLim)
            elseif j>i % Lower diagonal: scatter plot
                idx = 1:17:length(dataA);
                plot(dataA(idx,j), dataA(idx,i), 'o', 'markeredgecolor', Colors(1,:), 'markerfacecolor', Colors(1,:), 'markersize', 1)
                plot(dataB(idx,j), dataB(idx,i), 'o', 'markeredgecolor', Colors(4,:), 'markerfacecolor', Colors(4,:), 'markersize', 1)
            else % Upper diagonal: correlation coefficient
                correlCoefA = coefMatA(i,j);
                correlCoefB = coefMatB(i,j);
                text(0.5, 0.65, sprintf('No BC: %1.3f', correlCoefA), 'Color', Colors(1,:), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7)
                text(0.5, 0.45, sprintf('Yes BC: %1.3f', correlCoefB), 'Color', Colors(4,:), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7)
                xlim([0, 1])
                ylim([0, 1])
            end
            % text(0.5, 0.5, sprintf('n=%1.0f, i=%1.0f, j=%1.0f\n%s x %s', n, i, j, labels{j}, labels{i}))
            % Next subplot
            n = n + 1;
        end
    end 
    fName = ['Example2_FixedRecPlate/' plotName];
    print(gcf,fName,'-dpdf')
end

%%
function histo_plot(trueD,D_MLEa,histD_a,D_MLEb,histD_b,Plot_Name)
    Colors=[57  106 177
        218 124 48
        62  150 81
        204 37  41
        83  81  84
        107 76  154
        146 36  40
        148 139 61]./255;

    FigWidth=9;nby=1; nbx=1;spacey=1; spacex=0.7; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
    FigDepth=4+1.1;
    FontAxis=7; FontLabel=9; FontLegend=7;
    f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
    [ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);
    ax=axes('position',positions{1},'Layer','top'); hold on; box on; grid on;
    xPlot = 0.8:0.001:1.2; 
    numBins = 50;    
    h1 = histogram(histD_a./trueD, 'numbins', numBins, 'edgecolor', 'none', 'facecolor' ,Colors(1,:), 'facealpha', 0.3, 'Normalization','pdf');
    p1 = fitdist(histD_a./trueD, 'Normal'); p1 = pdf(p1, xPlot); p1 = plot(xPlot, p1, '-', 'color', Colors(1,:));   
    h2 = histogram(histD_b./trueD, 'numbins', numBins, 'edgecolor', 'none', 'facecolor' ,Colors(4,:), 'facealpha', 0.3, 'Normalization','pdf');
    p2 = fitdist(histD_b./trueD, 'Normal'); p2 = pdf(p2, xPlot); p2 = plot(xPlot, p2, '-', 'color', Colors(4,:));
     
    plotMaxY = ax.YLim(2);
    m1 = plot([D_MLEa./trueD, D_MLEa./trueD], [0, plotMaxY], '--', 'color', Colors(1,:));
    m2 = plot([D_MLEb./trueD, D_MLEb./trueD], [0, plotMaxY], '--', 'color', Colors(4,:));
   
    xlim(minmax(xPlot));
    ylim([0 plotMaxY]);
    xticks([0.8:0.1:1.2]);
    xtickformat('%.1f');  
    ytickformat('%.0f');
    set(gca,'FontSize',FontAxis);
    ylabel('$p(D)$ [-]','Interpreter', 'latex','FontSize',FontLabel)
    xlabel('$D/D_{\mathrm{true}}$ [-]','Interpreter', 'latex','FontSize',FontLabel)
    l=legend([p1 m1 p2 m2],{'No BC - MCMC','No BC - MLE','Yes BC - MCMC','Yes BC - MLE'},'location','best');
    set(l,'FontSize',FontLegend);
    fName = ['Example2_FixedRecPlate/' Plot_Name];
    print(gcf,fName,'-dpdf')
    hold off
end

%%
function cont_plot(x_ana,y_ana,f_ana,x_train,y_train,x_pred,y_pred,f_pred1,f_pred2,x_bc,y_bc,clim,ctik,cbar_label,Plot_Name)
Fact=0.7;
FigWidth=18.48*Fact;nby=1; nbx=3;spacey=1; spacex=0.7*Fact; leftmargin=1.2; rightmargin=0.5; topmargin=0.1; bottommargin=1.8;
FigDepth=(FigWidth-nbx*spacex)/3+bottommargin+topmargin;
FontAxis=7; FontLabel=9; FontLegend=7;NLvls=10;

f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

max_val_ana=max(max(abs(f_ana)));
ax=axes('position',positions{1},'Layer','top');
contourf(x_ana,y_ana,f_ana./max_val_ana,NLvls); hold on
colormap('coolwarm');
caxis(clim)
if ~isempty(x_train); plot(x_train,y_train,'xr'); end
if ~isempty(x_bc); plot(x_bc,y_bc,'ok'); end
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$y/b$ [-]','Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
box on; grid off; shading flat; axis equal

ax=axes('position',positions{2},'Layer','top');
contourf(x_pred,y_pred,f_pred1./max_val_ana,NLvls); hold on
colormap('coolwarm');
caxis(clim)
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('');
set(gca,'YTickLabel','')   
xtickformat('%.1f');
box on; grid off; shading flat; axis equal

ax=axes('position',positions{3},'Layer','top');
contourf(x_pred,y_pred,f_pred2./max_val_ana,NLvls); hold on
colormap('coolwarm');
caxis(clim)
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('');
set(gca,'YTickLabel','')   
xtickformat('%.1f');
box on; grid off; shading flat; axis equal

axes('Position', [positions{2,1}(1)+positions{2,1}(3)/2-0.246 (positions{1,1}(2)-0.325) positions{2,1}(3)*1.815 0.638], 'Visible', 'off');
cb=colorbar ('southoutside');
colormap('coolwarm');
set(gca,'Clim',clim, 'FontSize',FontAxis);set(cb,'XTick',ctik)
ylabel(cb,cbar_label, 'Interpreter', 'LaTeX','FontSize',FontLabel)
cb.Ruler.TickLabelFormat = '%.1f';

fName = ['Example2_FixedRecPlate/' Plot_Name];
print(gcf,fName,'-dpdf')
hold off
end

%%
function line_plot(varargin)

x_ana=varargin{1}; y_ana=varargin{2};f_ana=varargin{3};f_ana_p2=varargin{4};
x_train=varargin{5};y_train=varargin{6};f_train=varargin{7};f_train_p2=varargin{8};
x_pred=varargin{9};y_pred=varargin{10};f_pred1=varargin{11};f_pred1_std=varargin{12};f_pred2=varargin{13};f_pred2_std=varargin{14};
f_pred1_p2=varargin{15};f_pred1_std_p2=varargin{16};f_pred2_p2=varargin{17};f_pred2_std_p2=varargin{18};
EdgePoint=varargin{19};
BC=varargin{20};
xlab=varargin{21};ylab=varargin{22};xlab_p2=varargin{23};ylab_p2=varargin{24};LIM=varargin{25};TICK=varargin{26};
LegendNames=varargin{27};LegendNames1=varargin{28};
Plot_Name=varargin{29};

Colors=[57  106 177
        218 124 48
        62  150 81
        204 37  41
        83  81  84
        107 76  154
        146 36  40
        148 139 61]./255;

FigWidth=18.48;nby=1; nbx=2;spacey=1; spacex=2.18; leftmargin=1.2; rightmargin=0.5;topmargin=0.2; bottommargin=1.1;
FigDepth=4+1.3;
FontAxis=7; FontLabel=9; FontLegend=7;
  
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

if EdgePoint(1)==1
P=1; P_plot=1; P_pred=1;
elseif EdgePoint(1)==0.5
P=floor((length(y_train)-1)/4+1); 
P_plot=floor((length(y_ana)-1)/4+1);
P_pred=floor((length(y_pred)-1)/4+1);
else
P=floor((length(y_train)-1)/2+1); 
P_plot=floor((length(y_ana)-1)/2+1);
P_pred=floor((length(y_pred)-1)/2+1);
end
max_val_ana=max(max(abs(f_ana)));

f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
ax=axes('position',positions{1},'Layer','top');hold on;
h10=fill([x_pred; flip(x_pred)], [f_pred1(P_pred,:)'+2.56.*f_pred1_std(P_pred,:)'; flip(f_pred1(P_pred,:)'-2.56*f_pred1_std(P_pred,:)')]./max_val_ana,Colors(1,:),'LineStyle','none','facealpha',.3);
h11=fill([x_pred; flip(x_pred)], [f_pred2(P_pred,:)'+2.56.*f_pred2_std(P_pred,:)'; flip(f_pred2(P_pred,:)'-2.56*f_pred2_std(P_pred,:)')]./max_val_ana,Colors(4,:),'LineStyle','none','facealpha',.3);
h1=plot(x_ana,f_ana(P_plot,:)./max_val_ana,'-','color','k');
if ~isempty(f_train); h2=plot(x_train,f_train(P,:)./max_val_ana,'xr'); end
h3=plot(x_pred,f_pred1(P_pred,:)'./max_val_ana,'color',Colors(1,:));
h4=plot(x_pred,f_pred2(P_pred,:)'./max_val_ana,'color',Colors(4,:));
if ~isempty(BC)&&~isempty(BC{1}); h5=plot(BC{1}(1,:),BC{1}(2,:),'ok','MarkerSize',5);end
xlim(LIM{1});ylim(LIM{2});yticks(TICK{1});
set(gca, 'FontSize',FontAxis);
xlabel(xlab,'Interpreter', 'latex','FontSize',FontLabel)
ylabel(ylab,'Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
if ~isempty(f_train)&&~isempty(BC)&&~isempty(BC{1})
    l=legend([h1 h5 h2 h3 h4],LegendNames,'location','best');
elseif isempty(f_train)&&~isempty(BC)&&~isempty(BC{1})
    l=legend([h1 h5 h3 h4],LegendNames,'location','best');
elseif ~isempty(f_train)
    l=legend([h1 h2 h3 h4],LegendNames,'location','best');    
else
    l=legend([h1 h3 h4],LegendNames,'location','best');
end
set(l,'FontSize',FontLegend);
grid on; box on;

if EdgePoint(2)==1
P=1; P_plot=1; P_pred=1;
elseif EdgePoint(2)==0.5
P=floor((length(y_train)-1)/4+1); 
P_plot=floor((length(y_ana)-1)/4+1);
P_pred=floor((length(y_pred)-1)/4+1);
else
P=floor((length(y_train)-1)/2+1); 
P_plot=floor((length(y_ana)-1)/2+1);
P_pred=floor((length(y_pred)-1)/2+1);
end
max_val_ana=max(max(abs(f_ana_p2)));

ax=axes('position',positions{2},'Layer','top');hold on;
h10=fill([x_pred; flip(x_pred)], [f_pred1_p2(P_pred,:)'+2.56.*f_pred1_std_p2(P_pred,:)'; flip(f_pred1_p2(P_pred,:)'-2.56*f_pred1_std_p2(P_pred,:)')]./max_val_ana,Colors(1,:),'LineStyle','none','facealpha',.3);
h11=fill([x_pred; flip(x_pred)], [f_pred2_p2(P_pred,:)'+2.56.*f_pred2_std_p2(P_pred,:)'; flip(f_pred2_p2(P_pred,:)'-2.56*f_pred2_std_p2(P_pred,:)')]./max_val_ana,Colors(4,:),'LineStyle','none','facealpha',.3);
h1=plot(x_ana,f_ana_p2(P_plot,:)./max_val_ana,'-','color','k');
if ~isempty(f_train_p2); h2=plot(x_train,f_train_p2(P,:)./max_val_ana,'xr'); end
h3=plot(x_pred,f_pred1_p2(P_pred,:)'./max_val_ana,'color',Colors(1,:));
h4=plot(x_pred,f_pred2_p2(P_pred,:)'./max_val_ana,'color',Colors(4,:));
if ~isempty(BC)&&~isempty(BC{2}); h5=plot(BC{2}(1,:),BC{2}(2,:),'ok','MarkerSize',5);end
xlim(LIM{3});ylim(LIM{4});yticks(TICK{2});
set(gca, 'FontSize',FontAxis);
xlabel(xlab_p2,'Interpreter', 'latex','FontSize',FontLabel)
ylabel(ylab_p2,'Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
if ~isempty(f_train_p2)&&~isempty(BC)&&~isempty(BC{2})
    l=legend([h1 h5 h2 h3 h4],LegendNames1,'location','best');
elseif isempty(f_train_p2)&&~isempty(BC)&&~isempty(BC{2})
    l=legend([h1 h5 h3 h4],LegendNames1,'location','best');
elseif ~isempty(f_train_p2)
    l=legend([h1 h2 h3 h4],LegendNames1,'location','best');    
else
    l=legend([h1 h3 h4],LegendNames1,'location','best');
end
set(l,'FontSize',FontLegend);
grid on; box on;

fName = ['Example2_FixedRecPlate/' Plot_Name];
print(gcf,fName,'-dpdf')

end

%%
function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)
    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
    for i=1:nbx
       for j=1:nby
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           yfirst=bottommargin+(j-1.0)*(subysize+spacey);
           positions{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
       end
    end
end
