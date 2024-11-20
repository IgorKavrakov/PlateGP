function Example1b_Plots(D_Pred_Standard,D_Pred_MCMC,D)
    %%
    close all;

    %% Bar plot
    FontAxis=7; FontLabel=9; FontLegend=7;
    FigWidth=9; nbx=1;nby=1; leftmargin=1.2; rightmargin=0.5; topmargin=0.1; bottommargin=0.9; spacex=0;spacey=0;
    FigDepth=4+1;
    [ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

    Colors=[57  106 177
        204 37  41
        62  150 81
        218 124 48]./255;
    SNR_Loop=[5 10 20 100];
    
    FigType = 2;

        for i=1:size(D_Pred_Standard,3)
            for j=1:size(D_Pred_Standard,1)
                Rows=D_Pred_Standard(j,:,i)./D>0.2; %Remove optimiser stuck samples
                fprintf('MLE - Learn %d, SNR %d:\n',i,SNR_Loop(j))
                fprintf('Total of %f percent are disregarded when computing the mean of the Monte Carlo analysis\n \n', (length(find(Rows==0)))/1000*100)
                QuantL(j,1+2*(i-1)) = quantile(D_Pred_Standard(j,Rows,i)./D, 0.25);
                QuantH(j,1+2*(i-1)) = quantile(D_Pred_Standard(j,Rows,i)./D, 0.75);
                Mins(j,1+2*(i-1))=min(D_Pred_Standard(j,Rows,i))./D;
                Maxs(j,1+2*(i-1))=max(D_Pred_Standard(j,Rows,i))./D;
                Median(j,1+2*(i-1))=median(D_Pred_Standard(j,Rows,i))./D;
                Mean(j,1+2*(i-1))=mean(D_Pred_Standard(j,Rows,i))./D;
                STD(j,1+2*(i-1))=std(D_Pred_Standard(j,Rows,i))./D;
            end
        end

        for i=1:size(D_Pred_MCMC,3)
            for j=1:size(D_Pred_MCMC,1)
                QuantL(j,2*i) = quantile(D_Pred_MCMC(j,Rows,i)./D, 0.25);
                QuantH(j,2*i) = quantile(D_Pred_MCMC(j,Rows,i)./D, 0.75);
                Mins(j,2*i)=min(D_Pred_MCMC(j,Rows,i))./D;
                Maxs(j,2*i)=max(D_Pred_MCMC(j,Rows,i))./D;
                Median(j,2*i)=median(D_Pred_MCMC(j,Rows,i))./D;
                Mean(j,2*i)=mean(D_Pred_MCMC(j,Rows,i))./D;
                STD(j,2*i)=std(D_Pred_MCMC(j,Rows,i))./D;
            end
        end

        f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
        ax=axes('position',positions{1},'Layer','top');hold on;      
        ngroups = size(Mean, 1);
        nbars = size(Mean, 2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        barwidth = (groupwidth/nbars)*0.5*0.75;
        plot([0,5], [1,1], '-k');
        colList = [1,1,2,2,3,3]; colAlphaList = [1,0.4,1,0.4,1,0.4];
            b = [];
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                for j = 1:ngroups
                    thisX = [x(j)-barwidth, x(j)+barwidth, x(j)+barwidth, x(j)-barwidth];
                    thisY = [QuantL(j,i),QuantL(j,i),QuantH(j,i),QuantH(j,i)];
                    b(i) = fill(thisX, thisY, Colors(colList(i),:), 'edgecolor', 'none', 'facealpha', colAlphaList(i), 'displayName', sprintf('L%1.0f', i));
                    errorbar(x(j),Median(j,i),abs(Mins(j,i)-1),abs(Maxs(j,i)-1), '-','color',Colors(colList(i),:));
                    plot(x(j)+[-0.8, 0.8]*barwidth, Median(j,i)*[1,1], '-k')
                end
            end
        l=legend([b(1) b(3) b(5)],{'L1','L2','L3'},'location','southeast');
        set(l,'FontSize',FontLegend);
        xticks([1 2 3 4]);xticklabels([5 10 20 100]); xlim([0.5, 4.5])
        if FigType == 1
            ylim([0 1.6]); yticks([0:0.2:1.6]);
        else
            ylim([0.4 1.8]); yticks([0.4:0.2:1.8]);
        end
        set(gca,'FontSize',FontAxis);
        xlabel('Signal-to-noise (SNR) Ratio [-]')
        ylabel('$\overline{D}/D_{\mathrm{true}}$,         $\hat{D}/D_{\mathrm{true}}$ [-]','Interpreter', 'latex','FontSize',FontLabel)
        ytickformat('%.1f')
        box on;


    print(gcf,'Example1_SSRecPlate/Ex1b_D_ParStudySNR','-dpdf')

end

%% Plot function
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

