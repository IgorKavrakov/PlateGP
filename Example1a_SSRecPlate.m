%% READ-ME
% Stochastic Inference of Plate Bending from Heterogeneous Data: Physics-informed Gaussian Processes via Kirchhoff-Love theory
% 
% Cite: Kavrakov, I., Tondo. G. R., and Morgenthal, G. 2024. Stochastic Inference of Plate Bending from 
% Heterogeneous Data: Physics-informed Gaussian Processes via Kirchhoff-Love theory. 
% Journal of Engineering Mechanics, XXX, XXX, DOI:XXX.
% 
% PlateGP is a Matlab-based computer code for the inverse problem of obtaining plate stiffness values 
% from measurement data, and physics-informed forward solutions of plate physical responses. 
% The scripts Example1a_SSRecPlate.m and Example1b_SSRecPlate_SNR.m contain the 
% verification of the model for a rectangular simply supported plate, while Example2_FixedRecPlate.m implements a 
% plate model in fixed support conditions. The folder GP_Plate_Model includes an implementation of Gaussian Processes, 
% while the Example folders contain the details for the analytical flat plate analysis and plots. The details on of 
% the utilized methods are given in the abovementioned article.
% 
% The accepted manuscript can be found on arXiv: https://arxiv.org/abs/2405.12802
% 
% %%%%%%%%% COPYRIGHT NOTICE %%%%%%%%%
% 
% PlateGP is free software: you can redistribute it and/or modify it under the terms 
% of the GNU General Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% PlateGP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with PlateGP. 
% If not, see https://www.gnu.org/licenses/.
% 
% Copyright (c) Igor Kavrakov, Gledson Rodrigo Tondo, Guido Morgenthal 2024

%% Clear and add paths
clear all; close all; clc;matlabrc;
addpath(genpath('GP_Plate_Model'));     %Main dir for GP plate bending
addpath('Example1_SSRecPlate');         %Particular example dir: analytical solutions & plots

%% Input
rng(10); %Reproducability
%Analytical Plate properties
AP.nu=0.3;
AP.q0=1e3; %Sinusoidal load amplitude
AP.a=1; AP.b=1;   %Plate dim in x and y
AP.D=2.1*1e8*0.01^3/(12*(1-0.3^2)); %Thickness for learning!
AP.x=linspace(0,AP.a,5);AP.y=linspace(0,AP.b,5); %Learning coord
AP.x(1)=0.05*AP.a;AP.x(end)=0.95*AP.a;AP.y(1)=0.05*AP.b;AP.y(end)=0.95*AP.b; %Move outer coords
SNR=10; %Contaminate analytical data with specific SNR
AP.SNR.w=SNR; AP.SNR.Kx=SNR; AP.SNR.Ky=SNR; AP.SNR.Kxy=SNR; AP.SNR.p=SNR; 

%GP properties
GP.Par.Stabilizer=eps; %Jitter term.
GP.Par.noise=1; % Predict noise ?!
GP.Par.nu=AP.nu;

%% Optimizer properties
% Standard optimization
GP.Par.NEval=-10000; %Maximum function evaluations in training
% Bayesian model (Gledson Tondo)
GP.Par.chainLength=1e4; % MCMC chain length
GP.Par.burnIn = 2.5e3; % Discard initial MCMC steps, get stable part of the chain
GP.Par.thin = 2; % Discard intermediate steps, reduce dist. autocorrelation
GP.Par.jumpSize = 2e-2; % Jump size for parameter sampler, relative to initial parameter value

%% Plate analytical solution
[AP] = Example1_Analytical_SSRecPlate(AP);

%% Initiate GP Structure
Type='Init';
GP = GP_Plate_Pars(GP,Type); %Parser - Initialise data
GP(2)=GP(1); GP(3)=GP(1); %Save initial structures for 3 cases of learning

%% Displacement-based learning
%Learning
Type='Train'; %Training parser
Prop={'w','p'}; %Properties to get out of the Data structure (in this case, training data - AP)
GP(1)= GP_Plate_Pars(GP(1),Type,AP,Prop); % Parser - Feed in Data (AP).

hyp=log([1 1 1 exp(1) 10^-3 10^-3])';%Hyperparameters prior (D, A, l_x,l_y, 2 sigma_n);
fun=@(hyp)GP_Process_Opt(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,... % Objective function 
                         GP(1).Kernel,hyp,GP(1).Par.Stabilizer,GP(1).Par.noise);     
% Standard optimization
GP(1).Par.hyp= minimize(hyp,fun,GP(1).Par.NEval);%Minimize
GP(1).Par.D=exp(GP(1).Par.hyp(4));%Flexural stifness

% Bayesian optimization
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp); % log prior + log lokelihood
propRnd = @(hyp) normrnd(hyp, repmat(GP(1).Par.jumpSize, size(hyp))); % new parameter proposal
GP(1).Par.hyp_MCMC = metropolisHastings(hyp, GP(1).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(1).Par.burnIn, GP(1).Par.thin, true);
GP(1).Par.D_MCMC = mean(exp(GP(1).Par.hyp_MCMC(:,4)));%mean of flexural stifness
GP(1).Par.D_MCMC_stdev = std(exp(GP(1).Par.hyp_MCMC(:,4)));%standard deviation of flexural stifness
        

%Prediction
Type='Pred'; %Predition parser
Prop={'w','Kx','Ky','Kxy','Qx','Qy','Mx','My','Mxy'}; %Properties to predict
Pred.x=linspace(0,AP.a,21); 
Pred.y=linspace(0,AP.b,21); %Prediction (inference coord)
GP(1)= GP_Plate_Pars(GP(1),Type,Pred,Prop); % Parser - Feed in Prediction Data.

% Predictions: standard optimizer
[Pred.Mean,Pred.Var] = GP_Process(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,...
                                  GP(1).Pred.X,GP(1).Pred.Y,...
                                  GP(1).Kernel,GP(1).Par.hyp,GP(1).Par.Stabilizer,GP(1).Par.noise);
                              
% Predictions: Bayesian optimizer
[Pred.Mean_MCMC,Pred.Var_MCMC] = GP_Process(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,...
                                  GP(1).Pred.X,GP(1).Pred.Y,...
                                  GP(1).Kernel,GP(1).Par.hyp_MCMC,GP(1).Par.Stabilizer,GP(1).Par.noise);
                              
Type='Results'; %Predition parser
GP(1)= GP_Plate_Pars(GP(1),Type,Pred,Prop); % Parser - Feed in Results Data.

%% Curvature-based learning
Type='Train'; %Training parser
Prop={'Kx','Ky','Kxy','p'}; %Properties to get out of the Data structure (in this case, training data - AP)
GP(2)= GP_Plate_Pars(GP(2),Type,AP,Prop); % Parser - Feed in Data (AP).

hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3])';%Hyperparameters prior (D, A, l_x,l_y, 4 sigma_n);
fun=@(hyp)GP_Process_Opt(GP(2).Train.X,GP(2).Train.Y,GP(2).Train.f,... % Objective function 
                         GP(2).Kernel,hyp,GP(2).Par.Stabilizer,GP(2).Par.noise);    

% Standard optimization
GP(2).Par.hyp= minimize(hyp,fun,GP(2).Par.NEval);%Minimize
GP(2).Par.D=exp(GP(2).Par.hyp(4));%Flexural stifness

% Bayesian optimization
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp); % log prior + log lokelihood
propRnd = @(hyp) normrnd(hyp, repmat(GP(2).Par.jumpSize, size(hyp))); % new parameter proposal
GP(2).Par.hyp_MCMC = metropolisHastings(hyp, GP(2).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(2).Par.burnIn, GP(2).Par.thin, true);
GP(2).Par.D_MCMC = mean(exp(GP(2).Par.hyp_MCMC(:,4)));%mean of flexural stifness
GP(2).Par.D_MCMC_stdev = std(exp(GP(2).Par.hyp_MCMC(:,4)));%standard deviation of flexural stifness


%% Curvature-Displacement-based learning
%Learning
Type='Train'; %Training parser
Prop={'w','Kx','Ky','Kxy','p'}; %Properties to get out of the Data structure (in this case, training data - AP)
GP(3)= GP_Plate_Pars(GP(3),Type,AP,Prop); % Parser - Feed in Data (AP).

hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3 10^-3])';%Hyperparameters prior (D, A, l_x,l_y, 5 sigma_n);
fun=@(hyp)GP_Process_Opt(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,... % Objective function 
                         GP(3).Kernel,hyp,GP(3).Par.Stabilizer,GP(3).Par.noise);                     

% Standard optimization
GP(3).Par.hyp= minimize(hyp,fun,GP(3).Par.NEval);%Minimize
GP(3).Par.D=exp(GP(3).Par.hyp(4));%Flexural stifness

% Bayesian optimization
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp); % log prior + log lokelihood
propRnd = @(hyp) normrnd(hyp, repmat(GP(3).Par.jumpSize, size(hyp))); % new parameter proposal
GP(3).Par.hyp_MCMC = metropolisHastings(hyp, GP(3).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(3).Par.burnIn, GP(3).Par.thin, true);
GP(3).Par.D_MCMC = mean(exp(GP(3).Par.hyp_MCMC(:,4)));%mean of flexural stifness
GP(3).Par.D_MCMC_stdev = std(exp(GP(3).Par.hyp_MCMC(:,4)));%standard deviation of flexural stifness               


%Prediction
Type='Pred'; %Predition parser
Prop={'w','Kx','Ky','Kxy','Qx','Qy','Mx','My','Mxy'}; %Properties to predict
Pred.x=linspace(0,AP.a,21); 
Pred.y=linspace(0,AP.b,21); %Prediction (inference coord)
GP(3)= GP_Plate_Pars(GP(3),Type,Pred,Prop); % Parser - Feed in Prediction Data.

% Predictions: standard optimizer
[Pred.Mean,Pred.Var] = GP_Process(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,...
                                  GP(3).Pred.X,GP(3).Pred.Y,...
                                  GP(3).Kernel,GP(3).Par.hyp,GP(3).Par.Stabilizer,GP(3).Par.noise);
                              
% Predictions: Bayesian optimizer
[Pred.Mean_MCMC,Pred.Var_MCMC] = GP_Process(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,...
                                  GP(3).Pred.X,GP(3).Pred.Y,...
                                  GP(3).Kernel,GP(3).Par.hyp_MCMC,GP(3).Par.Stabilizer,GP(3).Par.noise);
                              
Type='Results'; %Predition parser
GP(3)= GP_Plate_Pars(GP(3),Type,Pred,Prop); % Parser - Feed in Results Data.
%% Save results
save('Example1_SSRecPlate/Example_1a.mat', 'GP','AP');

%% Plots
load('Example1_SSRecPlate/Example_1a.mat');
Example1a_Plots(GP,AP)
