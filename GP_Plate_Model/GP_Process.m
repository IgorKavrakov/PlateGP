function [Mean_Pred,Var,lik] = GP_Process(x,y,f,x_pred,y_pred,Kernel,hyp,Stabilizer,noise)
% Implements GP model

% By Igor Kavrakov

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of PlateGP.
%  PlateGP is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  PlateGP is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with PlateGP.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Gledson Rodrigo Tondo, Guido Morgenthal 2024
%This is implementation of GP regression by Igor Kavrakov based on 
%Rasmussen&Willams Gaussian Processes for Machine Learning (2006) (Algorithm 2.1)
%Adapted by Gledson Tondo for prediction with probabilistic parameters
isProbFlag = all(size(hyp) > 1);
if isProbFlag
    % Verbose
    % Predictions with probabilistic params
    % Allocate storage for results
    predTargets = fieldnames(x_pred); 
    numTargets = 0;
    for i = 1:length(predTargets) 
        if ~isempty(x_pred.(predTargets{i})) 
            numTargets = numTargets + 1; 
            numPredPoints = length(x_pred.(predTargets{i}));
        end 
    end
    numMC = 100; % Number of empirical MC iterations
    if numMC > size(hyp,1); numMC = size(hyp,1); end
    Mean_Pred = zeros(numTargets*numPredPoints, numMC);
    Var = zeros(numTargets*numPredPoints, numTargets*numPredPoints);
    lik = zeros(numMC, 1);
    % Create hyperparameter samples
    % Uniformly sampling from MCMC chain, we can also model params individually with some 
    % distribution and then sample from it, but this might imply losing some correlation
    idx = randi(size(hyp,1), [numMC, 1]);
    fprintf(1,'Bayesian prediction: %3d%%\n', 0)
    for i = 1:numMC
        [Mean_Pred(:,i),thisVar,lik(i)] = standardGPPreds(x,y,f,x_pred,y_pred,Kernel,hyp(idx(i),:)',Stabilizer,noise);
        Var = Var + thisVar;
        fprintf(1,'\b\b\b\b\b%3.0f%%\n', 100*i/numMC);
    end
    % Empirical prob predictions
    empVar = zeros(size(Var));
    tmpMean = mean(Mean_Pred, 2);
    for i = 1:numMC
        deltaMean = Mean_Pred(:,i) - tmpMean;
        empVar = empVar + deltaMean*deltaMean';
    end
    Mean_Pred = tmpMean;
    Var = Var./numMC + empVar./numMC;
else
    [Mean_Pred,Var,lik] = standardGPPreds(x,y,f,x_pred,y_pred,Kernel,hyp,Stabilizer,noise);
end
end

function [Mean_Pred,Var,lik] = standardGPPreds(x,y,f,x_pred,y_pred,Kernel,hyp,Stabilizer,noise)
    fy=[f.w;f.Rx;f.Ry;f.Kx;f.Ky;f.Kxy;f.p;f.Qx;f.Qy;f.Mx;f.My;f.Mxy]; %Reshape the output;
    Kern=Kern_Ensemble(Kernel.Cov,Kernel.BC,Kernel.nu,x,x,y,y,hyp,noise);%Get the kernel at x points
    Kern_predict=Kern_Ensemble(Kernel.Cov,[],Kernel.nu,x_pred,x_pred,y_pred,y_pred,hyp,0); %Get the kernel at prediction points
    Kern_predict_r=Kern_Ensemble(Kernel.Cov,[],Kernel.nu,x_pred,x,y_pred,y,hyp,0); %Get kernel at cross points

    Kern=Kern+eye(size(Kern))*Stabilizer; %Stabilzer (i.e. jitter term)

    L_Kern=chol(Kern,'lower'); %Compute cholesky (for inversion - faster). Check Rasmussen
    alpha=L_Kern'\(L_Kern\fy);       %Predictive mean
    Mean_Pred=Kern_predict_r*alpha; %Predictive mean
    v=L_Kern\Kern_predict_r';       %Predictive varianse
    Var=Kern_predict-v'*v;          %Predictive variances
    lik=-1/2*fy'*alpha-sum(diag(L_Kern))-length(fy)/2*log(2*pi); %Likelihood
end

