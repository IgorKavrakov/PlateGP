%% READ-ME
% Stochastic Inference of Plate Bending from Heterogeneous Data: Physics-informed Gaussian Processes via Kirchhoff-Love theory
% 
% Cite: Kavrakov, I., Tondo. G. R., and Morgenthal, G. 2025. Stochastic Inference of Plate Bending from 
% Heterogeneous Data: Physics-informed Gaussian Processes via Kirchhoff-Love theory. 
% Journal of Engineering Mechanics 151 (4), DOI: 10.1061/JENMDT.EMENG-7558.
% 
% PlateGP is a Matlab-based computer code for the inverse problem of obtaining plate stiffness values 
% from measurement data, and physics-informed forward solutions of plate physical responses. 
% The scripts Example1a_SSRecPlate.m and Example1b_SSRecPlate_SNR.m contain the 
% verification of the model for a rectangular simply supported plate, while Example2_FixedRecPlate.m implements a 
% plate model in fixed support conditions. The folder GP_Plate_Model includes an implementation of Gaussian Processes, 
% while the Example folders contain the details for the analytical flat plate analysis and plots. The details on of 
% the utilized methods are given in the abovementioned article.
% 
% The accepted manuscript can be found at https://doi.org/10.1061/JENMDT.EMENG-7558. 
% A preprint is available at arXiv: https://arxiv.org/abs/2405.12802
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
% Copyright (c) Igor Kavrakov, Gledson Rodrigo Tondo, Guido Morgenthal 2025

%% Clear and add paths
clear; close all; clc; matlabrc;

% Add paths to model and example
addpath(genpath('GP_Plate_Model'));
addpath('Example1_SSRecPlate');

% For reproducibility
rng(10);

%% Analytical plate properties and response
AP.nu=0.3;                                          % Poisson ratio
AP.q0=1e3;                                          % Load amplitude
AP.a=1; AP.b=1;                                     % Plate dimensions in x and y
AP.D=2.1*1e8*0.01^3/(12*(1-0.3^2));                 % Thickness for learning
AP.x=linspace(0,AP.a,5); AP.y=linspace(0,AP.b,5);   % Coordinates for training examples
AP.x(1)=0.05*AP.a;AP.x(end)=0.95*AP.a;              % Move outer coordinates
AP.y(1)=0.05*AP.b;AP.y(end)=0.95*AP.b; 
SNR=10;                                             % Contaminate analytical data with specific SNR
AP.SNR.w=SNR; AP.SNR.Kx=SNR; AP.SNR.Ky=SNR; 
AP.SNR.Kxy=SNR; AP.SNR.p=SNR; 

% Compute analytical solution
[AP] = Example1_Analytical_SSRecPlate(AP);

%% Gaussian process properties
GP.Par.Stabilizer=eps;      % Jitter term.
GP.Par.noise=1;             % Predict noise
GP.Par.nu=AP.nu;            % Plate Poisson ratio

%% Optimizer properties
% Standard optimization: maximum function evaluations in training
GP.Par.NEval=-10000;

% Bayesian optimization via MCMC
GP.Par.chainLength=1e4; % MCMC chain length
GP.Par.burnIn = 2.5e3;  % Discard initial MCMC steps, get stable part of the chain
GP.Par.thin = 2;        % Discard intermediate steps, reduce dist. autocorrelation
GP.Par.jumpSize = 2e-2; % Jump size for parameter sampler, relative to initial parameter value

%% Initiate GP Structure
Type='Init';
GP = GP_Plate_Pars(GP,Type);    % Parser - Create empty GP framework
GP(2)=GP(1); GP(3)=GP(1);       % Copy base model for 3 learning examples

%% Example 1: Displacement-based learning
% Add training data to the GP model
Type='Train';       % Training parser
Prop={'w','p'};     % Properties to add as training data
GP(1)= GP_Plate_Pars(GP(1),Type,AP,Prop);

% Initial value of hyperparameters (D, A, lx, ly, noise)
hyp=log([1 1 1 exp(1) 10^-3 10^-3])';

% Define training objective function
fun=@(hyp)GP_Process_Opt(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,...
                         GP(1).Kernel,hyp,GP(1).Par.Stabilizer,GP(1).Par.noise);     

% Standard optimization: gradient descent
GP(1).Par.hyp= minimize(hyp,fun,GP(1).Par.NEval);

% Get value of plate stiffness
GP(1).Par.D=exp(GP(1).Par.hyp(4));

% Bayesian optimization
% Prior model
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp);
% Model for proposal of new parameters
propRnd = @(hyp) normrnd(hyp, repmat(GP(1).Par.jumpSize, size(hyp)));

% Run Bayesian optimisation
GP(1).Par.hyp_MCMC = metropolisHastings(hyp, GP(1).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(1).Par.burnIn, GP(1).Par.thin, true);
                               
% Get mean and standard deviation of plate stiffness                               
GP(1).Par.D_MCMC = mean(exp(GP(1).Par.hyp_MCMC(:,4)));
GP(1).Par.D_MCMC_stdev = std(exp(GP(1).Par.hyp_MCMC(:,4)));
        

% Prediction
% Define properties to predict
Type='Pred'; 
Prop={'w','Kx','Ky','Kxy','Qx','Qy','Mx','My','Mxy'}; 

% Define locations to predict
Pred.x=linspace(0,AP.a,21); 
Pred.y=linspace(0,AP.b,21);

% Parser - Feed in prediction data
GP(1)= GP_Plate_Pars(GP(1),Type,Pred,Prop); 

% Predictions: standard optimizer
[Pred.Mean,Pred.Var] = GP_Process(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,...
                                  GP(1).Pred.X,GP(1).Pred.Y,...
                                  GP(1).Kernel,GP(1).Par.hyp,GP(1).Par.Stabilizer,GP(1).Par.noise);
                              
% Predictions: Bayesian optimizer
[Pred.Mean_MCMC,Pred.Var_MCMC] = GP_Process(GP(1).Train.X,GP(1).Train.Y,GP(1).Train.f,...
                                  GP(1).Pred.X,GP(1).Pred.Y,...
                                  GP(1).Kernel,GP(1).Par.hyp_MCMC,GP(1).Par.Stabilizer,GP(1).Par.noise);

% Parser: store results
Type='Results';
GP(1)= GP_Plate_Pars(GP(1),Type,Pred,Prop);

%% Example 2: Curvature-based learning
% Add training data to the GP model
Type='Train';               %Training parser
Prop={'Kx','Ky','Kxy','p'}; % Properties to add as training data
GP(2)= GP_Plate_Pars(GP(2),Type,AP,Prop);

% Initial value of hyperparameters (D, A, lx, ly, noise)
hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3])';

% Define training objective function
fun=@(hyp)GP_Process_Opt(GP(2).Train.X,GP(2).Train.Y,GP(2).Train.f,... % Objective function 
                         GP(2).Kernel,hyp,GP(2).Par.Stabilizer,GP(2).Par.noise);    

% Standard optimization
GP(2).Par.hyp= minimize(hyp,fun,GP(2).Par.NEval);

% Get value of plate stiffness
GP(2).Par.D=exp(GP(2).Par.hyp(4));

% Bayesian optimization
% Prior model
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp);
% Model for proposal of new parameters
propRnd = @(hyp) normrnd(hyp, repmat(GP(2).Par.jumpSize, size(hyp)));

% Run Bayesian optimisation
GP(2).Par.hyp_MCMC = metropolisHastings(hyp, GP(2).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(2).Par.burnIn, GP(2).Par.thin, true);
                               
% Get mean and standard deviation of plate stiffness   
GP(2).Par.D_MCMC = mean(exp(GP(2).Par.hyp_MCMC(:,4)));
GP(2).Par.D_MCMC_stdev = std(exp(GP(2).Par.hyp_MCMC(:,4)));

%% Example 3: Curvature-Displacement-based learning
% Add training data to the GP model
Type='Train';                   % Training parser
Prop={'w','Kx','Ky','Kxy','p'}; % Properties to get out of the Data structure
GP(3)= GP_Plate_Pars(GP(3),Type,AP,Prop);

% Initial value of hyperparameters (D, A, lx, ly, noise)
hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3 10^-3])';

% Define training objective function
fun=@(hyp)GP_Process_Opt(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,... % Objective function 
                         GP(3).Kernel,hyp,GP(3).Par.Stabilizer,GP(3).Par.noise);                     

% Standard optimization
GP(3).Par.hyp= minimize(hyp,fun,GP(3).Par.NEval);

% Get value of plate stiffness
GP(3).Par.D=exp(GP(3).Par.hyp(4));

% Bayesian optimization
% Prior model
funMCMC = @(hyp) GP_HypsPrior(hyp) - fun(hyp);
% Model for proposal of new parameters
propRnd = @(hyp) normrnd(hyp, repmat(GP(3).Par.jumpSize, size(hyp)));

% Run Bayesian optimisation
GP(3).Par.hyp_MCMC = metropolisHastings(hyp, GP(3).Par.chainLength, funMCMC, propRnd, [], ....
                                   GP(3).Par.burnIn, GP(3).Par.thin, true);
                               
% Get mean and standard deviation of plate stiffness                              
GP(3).Par.D_MCMC = mean(exp(GP(3).Par.hyp_MCMC(:,4)));
GP(3).Par.D_MCMC_stdev = std(exp(GP(3).Par.hyp_MCMC(:,4)));            


% Prediction
% Define properties to predict
Type='Pred'; 
Prop={'w','Kx','Ky','Kxy','Qx','Qy','Mx','My','Mxy'};

% Define locations to predict
Pred.x=linspace(0,AP.a,21); 
Pred.y=linspace(0,AP.b,21);

% Parser - Feed in prediction data
GP(3)= GP_Plate_Pars(GP(3),Type,Pred,Prop); 

% Predictions: standard optimizer
[Pred.Mean,Pred.Var] = GP_Process(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,...
                                  GP(3).Pred.X,GP(3).Pred.Y,...
                                  GP(3).Kernel,GP(3).Par.hyp,GP(3).Par.Stabilizer,GP(3).Par.noise);
                              
% Predictions: Bayesian optimizer
[Pred.Mean_MCMC,Pred.Var_MCMC] = GP_Process(GP(3).Train.X,GP(3).Train.Y,GP(3).Train.f,...
                                  GP(3).Pred.X,GP(3).Pred.Y,...
                                  GP(3).Kernel,GP(3).Par.hyp_MCMC,GP(3).Par.Stabilizer,GP(3).Par.noise);

% Parser: store results
Type='Results';
GP(3)= GP_Plate_Pars(GP(3),Type,Pred,Prop);

%% Save results to file
save('Example1_SSRecPlate/Example_1a.mat', 'GP','AP');

%% Plots
load('Example1_SSRecPlate/Example_1a.mat');
Example1a_Plots(GP,AP)
