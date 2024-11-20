function Example1a_Plots(GP,AP)
%% This is a specific function for plots in the paper
close all; clc;

%% Learned flexural rigidity
fprintf('Flexural Rigidity:\n')
fprintf('Learn 1: (Standard) D/D_true= %2.3f \t\t (MCMC) D/D_true= %2.3f\n', GP(1).Par.D/AP.D, GP(1).Par.D_MCMC/AP.D)
fprintf('Learn 2: (Standard) D/D_true= %2.3f \t\t (MCMC) D/D_true= %2.3f\n', GP(2).Par.D/AP.D, GP(2).Par.D_MCMC/AP.D)
fprintf('Learn 3: (Standard) D/D_true= %2.3f \t\t (MCMC) D/D_true= %2.3f\n', GP(3).Par.D/AP.D, GP(2).Par.D_MCMC/AP.D)

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
    [AP] = Example1_Analytical_SSRecPlate(AP);

    %% Figure 1 - Contour Displacements w
    cont_plot(AP.x_mesh,AP.y_mesh,AP.w,...
              GP(1).Train.x_mesh.w,GP(1).Train.y_mesh.w,...
              GP(1).Pred.x_mesh.w,GP(1).Pred.y_mesh.w,GP(1).Pred.(predMeshMean).w,GP(3).Pred.(predMeshMean).w,...
              '$w/\max |w_{\mathrm{true}}|$',sprintf('Ex1a_Contour_w%s', subName))

    %% Figure 2 - Contour Curvature Kx
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Kx,...
              GP(3).Train.x_mesh.Kx,GP(3).Train.y_mesh.Kx,...
              GP(1).Pred.x_mesh.Kx,GP(1).Pred.y_mesh.Kx,GP(1).Pred.(predMeshMean).Kx,GP(3).Pred.(predMeshMean).Kx,...
              '$\kappa_x/\max|\kappa_{x,\mathrm{true}}|$',sprintf('Ex1a_Contour_Kx%s', subName))

    %% Figure 3 - Contour Curvature Ky
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Ky,...
              GP(3).Train.x_mesh.Ky,GP(3).Train.y_mesh.Ky,...
              GP(1).Pred.x_mesh.Ky,GP(1).Pred.y_mesh.Ky,GP(1).Pred.(predMeshMean).Ky,GP(3).Pred.(predMeshMean).Ky,...
              '$\kappa_y/\max|\kappa_{y,\mathrm{true}}|$',sprintf('Ex1a_Contour_Ky%s', subName))

    %% Figure 4 - Contour Curvature Kxy
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Kxy,...
              GP(3).Train.x_mesh.Kxy,GP(3).Train.y_mesh.Kxy,...
              GP(1).Pred.x_mesh.Kxy,GP(1).Pred.y_mesh.Kxy,GP(1).Pred.(predMeshMean).Kxy,GP(3).Pred.(predMeshMean).Kxy,...
              '$\kappa_{xy}/\max|\kappa_{xy,\mathrm{true}}|$',sprintf('Ex1a_Contour_Kxy%s', subName))

    %% Figure 5 - Contour Shear Qx
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Qx,...
              [],[],...
              GP(1).Pred.x_mesh.Qx,GP(1).Pred.y_mesh.Qx,GP(1).Pred.(predMeshMean).Qx,GP(3).Pred.(predMeshMean).Qx,...
              '$Q_x/\max|Q_{x,\mathrm{true}}|$',sprintf('Ex1a_Contour_Qx%s', subName))      

    %% Figure 6 - Contour Shear Qy
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Qy,...
              [],[],...
              GP(1).Pred.x_mesh.Qy,GP(1).Pred.y_mesh.Qy,GP(1).Pred.(predMeshMean).Qy,GP(3).Pred.(predMeshMean).Qy,...
              '$Q_y/\max|Q_{y,\mathrm{true}}|$',sprintf('Ex1a_Contour_Qy%s', subName))          

    %% Figure 7 - Contour Moment Mx
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Mx,...
              [],[],...
              GP(1).Pred.x_mesh.Mx,GP(1).Pred.y_mesh.Mx,GP(1).Pred.(predMeshMean).Mx,GP(3).Pred.(predMeshMean).Mx,...
              '$M_x/\max|M_{x,\mathrm{true}}|$',sprintf('Ex1a_Contour_Mx%s', subName))         

    %% Figure 8 - Contour Moment My
    cont_plot(AP.x_mesh,AP.y_mesh,AP.My,...
              [],[],...
              GP(1).Pred.x_mesh.My,GP(1).Pred.y_mesh.My,GP(1).Pred.(predMeshMean).My,GP(3).Pred.(predMeshMean).My,...
              '$M_y/\max|M_{y,\mathrm{true}}|$',sprintf('Ex1a_Contour_My%s', subName))

    %% Figure 9 - Contour Moment Mxy
    cont_plot(AP.x_mesh,AP.y_mesh,AP.Mxy,...
              [],[],...
              GP(1).Pred.x_mesh.Mxy,GP(1).Pred.y_mesh.Mxy,GP(1).Pred.(predMeshMean).Mxy,GP(3).Pred.(predMeshMean).Mxy,...
              '$M_{xy}/\max|M_{xy,\mathrm{true}}|$',sprintf('Ex1a_Contour_Mxy%s', subName))    

    %% Figure 10 - Line plot Displacement & Curvature Kx
    line_plot(AP.x,AP.y,AP.w,AP.Kx,...
              GP(3).Train.x.w,GP(3).Train.y.w,GP(3).Train.f_mesh.w,GP(3).Train.f_mesh.Kx,...
              GP(3).Pred.x.w,GP(3).Pred.y.w,GP(1).Pred.(predMeshMean).w,GP(1).Pred.(predMeshStd).w,GP(3).Pred.(predMeshMean).w,GP(3).Pred.(predMeshStd).w,...
              GP(1).Pred.(predMeshMean).Kx,GP(1).Pred.(predMeshStd).Kx,GP(3).Pred.(predMeshMean).Kx,GP(3).Pred.(predMeshStd).Kx,...
              [0 0],...
              '$x/a$ [-]','$w/\max |w_{\mathrm{true}}|$ [-]','$x/a$ [-]','$\kappa_x/\max |\kappa_{x,\mathrm{true}}|$ [-]',...
              {[0 1],[-0.2 1.201],[0 1],[-0.2 1.201]},{-0.2:0.2:1.2,-0.2:0.2:1.2},...
              {'True','Observation L1/L3','Prediction L1','Prediction L3'},...
              {'True','Observation L3','Prediction L1','Prediction L3'},...          
              sprintf('Ex1a_CenterLine_wKx%s', subName));

    %% Figure 11 - Line plot Shear Qx & Moment Mx
    line_plot(AP.x,AP.y,AP.Qx,AP.Mx,...
              [],[],[],[],...
              GP(3).Pred.x.Qx,GP(3).Pred.y.Qx,GP(1).Pred.(predMeshMean).Qx,GP(1).Pred.(predMeshStd).Qx,GP(3).Pred.(predMeshMean).Qx,GP(3).Pred.(predMeshStd).Qx,...
              GP(1).Pred.(predMeshMean).Mx,GP(1).Pred.(predMeshStd).Mx,GP(3).Pred.(predMeshMean).Mx,GP(3).Pred.(predMeshStd).Mx,...
              [0 0],...
              '$x/a$ [-]','$Q_x/\max |Q_{x,\mathrm{true}}|$ [-]','$x/a$ [-]','$M_x/\max |M_{x,\mathrm{true}}|$ [-]',...
              {[0 1],[-1.201 1.201],[0 1],[-0.2 1.201]},{-1.2:0.4:1.2,-0.2:0.2:1.2},...
              {'True','Prediction L1','Prediction L3'},...
              {'True','Prediction L1','Prediction L3'},...         
              sprintf('Ex1a_CenterLine_QxMx%s', subName));      

    %% Figure 12 - Line plot Curvature Kxy & Moment Mxy
    line_plot(AP.x,AP.y,AP.Kxy,AP.Mxy,...
              GP(3).Train.x.Kxy,GP(3).Train.y.Kxy,GP(3).Train.f_mesh.Kxy,[],...
              GP(3).Pred.x.Kxy,GP(3).Pred.y.Kxy,GP(1).Pred.(predMeshMean).Kxy,GP(1).Pred.(predMeshStd).Kxy,GP(3).Pred.(predMeshMean).Kxy,GP(3).Pred.(predMeshStd).Kxy,...
              GP(1).Pred.(predMeshMean).Mxy,GP(1).Pred.(predMeshStd).Mxy,GP(3).Pred.(predMeshMean).Mxy,GP(3).Pred.(predMeshStd).Mxy,...
              [1 1],...
              '$x/a$ [-]','$\kappa_{xy}/\max |\kappa_{xy,\mathrm{true}}|$ [-]','$x/a$ [-]','$M_{xy}/\max |M_{xy,\mathrm{true}}|$ [-]',...
              {[0 1],[-1.21 1.21],[0 1],[-1.21 1.21]},{-1.2:0.4:1.2,-1.2:0.4:1.2},...
              {'True','Observation L3','Prediction L1','Prediction L3'},...
              {'True','Prediction L1','Prediction L3'},...         
              sprintf('Ex1a_CenterLine_KxyMxy%s', subName));      

end

    %% Figure 13 - Histogram of D for each learning type - only if MCMC
        histo_plot(AP.D,GP(1).Par.D,exp(GP(1).Par.hyp_MCMC(:,4)),...
                        GP(2).Par.D,exp(GP(2).Par.hyp_MCMC(:,4)),...
                        GP(3).Par.D,exp(GP(3).Par.hyp_MCMC(:,4)),...
                   sprintf('Ex1a_Histogram_D%s',subName))
end

%%
function histo_plot(trueD,D_MLEa,histD_a,D_MLEb,histD_b,D_MLEc,histD_c,Plot_Name)
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
    ax=axes('position',positions{1},'Layer','top'); hold on; box on;; grid on;
    xPlot = 0.8:0.001:1.2; 
    numBins = 50;    
    h1 = histogram(histD_a./trueD, 'numbins', numBins, 'edgecolor', 'none', 'facecolor' ,Colors(1,:), 'facealpha', 0.3, 'Normalization','pdf');
    p1 = fitdist(histD_a./trueD, 'Normal'); p1 = pdf(p1, xPlot); p1 = plot(xPlot, p1, '-', 'color', Colors(1,:));   
    h2 = histogram(histD_b./trueD, 'numbins', numBins, 'edgecolor', 'none', 'facecolor' ,Colors(3,:), 'facealpha', 0.3, 'Normalization','pdf');
    p2 = fitdist(histD_b./trueD, 'Normal'); p2 = pdf(p2, xPlot); p2 = plot(xPlot, p2, '-', 'color', Colors(3,:));
    h3 = histogram(histD_c./trueD, 'numbins', numBins, 'edgecolor', 'none', 'facecolor' ,Colors(4,:), 'facealpha', 0.3, 'Normalization','pdf');
    p3 = fitdist(histD_c./trueD, 'Normal'); p3 = pdf(p3, xPlot); p3 = plot(xPlot, p3, '-', 'color', Colors(4,:));
    
    plotMaxY = ax.YLim(2);
    m1 = plot([D_MLEa./trueD, D_MLEa./trueD], [0, plotMaxY], '--', 'color', Colors(1,:));
    m2 = plot([D_MLEb./trueD, D_MLEb./trueD], [0, plotMaxY], '--', 'color', Colors(3,:));
    m3 = plot([D_MLEc./trueD, D_MLEc./trueD], [0, plotMaxY], '--', 'color', Colors(4,:));
   
    xlim(minmax(xPlot));
    ylim([0 plotMaxY]);
    xticks([0.8:0.1:1.2]);
    xtickformat('%.1f');  
    ytickformat('%.0f');
    set(gca,'FontSize',FontAxis);
    ylabel('$p(D)$ [-]','Interpreter', 'latex','FontSize',FontLabel)
    xlabel('$D/D_{\mathrm{true}}$ [-]','Interpreter', 'latex','FontSize',FontLabel)
    l=legend([p1 m1 p2 m2 p3 m3],{'L1 - MCMC','L1 - MLE','L2 - MCMC','L2 - MLE','L3 - MCMC','L3 - MLE'},'location','best');
    set(l,'FontSize',FontLegend);
    fName = ['Example1_SSRecPlate/' Plot_Name];
    print(gcf,fName,'-dpdf')
    hold off
end

%%
function cont_plot(x_ana,y_ana,f_ana,x_train,y_train,x_pred,y_pred,f_pred1,f_pred2,cbar_label,Plot_Name)
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
if min(min(f_ana))>-0.01 caxis([-0.2 1.2]); else; caxis([-1.2 1.2]); end
if ~isempty(x_train); plot(x_train,y_train,'xr'); end
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$y/b$ [-]','Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
box on; grid off; shading flat; axis equal

ax=axes('position',positions{2},'Layer','top');
contourf(x_pred,y_pred,f_pred1./max_val_ana,NLvls); hold on
colormap('coolwarm');
if min(min(f_ana))>-0.01 caxis([-0.2 1.2]); else; caxis([-1.2 1.2]); end
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('');
set(gca,'YTickLabel','')   
xtickformat('%.1f');
box on; grid off; shading flat; axis equal

ax=axes('position',positions{3},'Layer','top');
contourf(x_pred,y_pred,f_pred2./max_val_ana,NLvls); hold on
colormap('coolwarm');
if min(min(f_ana))>-0.01 caxis([-0.2 1.2]); else; caxis([-1.2 1.2]); end
set(gca, 'FontSize',FontAxis);
xticks(0:0.2:1);xlabel('$x/a$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('');
set(gca,'YTickLabel','')   
xtickformat('%.1f');
box on; grid off; shading flat; axis equal

axes('Position', [positions{2,1}(1)+positions{2,1}(3)/2-0.246 (positions{1,1}(2)-0.325) positions{2,1}(3)*1.815 0.638], 'Visible', 'off');
cb=colorbar ('southoutside');
if min(min(f_ana))>-0.01 
set(gca,'Clim',[-0.2 1.2], 'FontSize',FontAxis);set(cb,'XTick',-0.2:0.2:1.2)
else
set(gca,'Clim',[-1.2 1.2], 'FontSize',FontAxis);set(cb,'XTick',-1.2:0.4:1.2)
end
ylabel(cb,cbar_label, 'Interpreter', 'LaTeX','FontSize',FontLabel)
cb.Ruler.TickLabelFormat = '%.1f';

fName = ['Example1_SSRecPlate/' Plot_Name];
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
xlab=varargin{20};ylab=varargin{21};xlab_p2=varargin{22};ylab_p2=varargin{23};LIM=varargin{24};TICK=varargin{25};
LegendNames=varargin{26};LegendNames1=varargin{27};
Plot_Name=varargin{28};

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

if EdgePoint(1)
P=1; P_plot=1; P_pred=1;
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
xlim(LIM{1});ylim(LIM{2});yticks(TICK{1});
set(gca, 'FontSize',FontAxis);
xlabel(xlab,'Interpreter', 'latex','FontSize',FontLabel)
ylabel(ylab,'Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
if ~isempty(f_train)
    l=legend([h1 h2 h3 h4],LegendNames,'location','best');
else
    l=legend([h1 h3 h4],LegendNames,'location','best');
end
set(l,'FontSize',FontLegend);
grid on; box on;

if EdgePoint(2)
P=1; P_plot=1; P_pred=1;
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
xlim(LIM{3});ylim(LIM{4});yticks(TICK{2});
set(gca, 'FontSize',FontAxis);
xlabel(xlab_p2,'Interpreter', 'latex','FontSize',FontLabel)
ylabel(ylab_p2,'Interpreter', 'latex','FontSize',FontLabel)
xtickformat('%.1f');ytickformat('%.1f');
if ~isempty(f_train_p2)
    l=legend([h1 h2 h3 h4],LegendNames1,'location','best');
else
    l=legend([h1 h3 h4],LegendNames1,'location','best');
end
set(l,'FontSize',FontLegend);
grid on; box on;

fName = ['Example1_SSRecPlate/' Plot_Name];
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
