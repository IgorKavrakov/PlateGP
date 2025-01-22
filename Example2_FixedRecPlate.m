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
addpath('Example2_FixedRecPlate');  

% For reproducibility
rng(10);

%% Analytical plate properties and response
AP.nu=0.3;                                          % Poisson ratio                 
AP.q0=1e3;                                          % Load amplitude
AP.a=1; AP.b=1;                                     % Plate dimensions in x and y
AP.D=2.1*1e8*0.01^3/(12*(1-0.3^2));                 % Thickness for learning
AP.x=linspace(0,AP.a,5);AP.y=linspace(0,AP.b,5);    % Coordinates for training examples
AP.x(1)=0.05*AP.a;AP.x(end)=0.95*AP.a;              % Move outer coordinates
AP.y(1)=0.05*AP.b;AP.y(end)=0.95*AP.b; 
SNR=10;                                             % Contaminate analytical data with specific SNR
AP.SNR.w=SNR; AP.SNR.Kx=SNR; AP.SNR.Ky=SNR; 
AP.SNR.Kxy=SNR; AP.SNR.p=SNR; 

% Compute analytical solution
[AP] = Example2_Analytical_FixedRecPlate(AP);

%% Gaussian process properties
GP.Par.Stabilizer=1e-6;         % Jitter term
GP.Par.noise=1;                 % Predict noise
GP.Par.nu=AP.nu;                % Plate Poisson ratio

%% Optimizer properties
% Standard optimization: maximum function evaluations in training
GP.Par.NEval=-200; 

% Bayesian optimization via MCMC
GP.Par.chainLength=5e4;     % MCMC chain length
GP.Par.burnIn = 1e4;        % Discard initial MCMC steps, get stable part of the chain
GP.Par.thin = 2;            % Discard intermediate steps, reduce dist. autocorrelation
GP.Par.jumpSize = 2e-2;     % Jump size for parameter sampler, relative to initial parameter value

%% Initiate GP Structure
Type='Init';
GP = GP_Plate_Pars(GP,Type);    % Parser - Create empty GP framework
GP(2)=GP(1);                    % Copy base model for 3 learning examples

%% Example 1: Curvature-Displacement-based learning - No Boundary condtions applied
% Add training data to the GP model
Type='Train';                       % Training parser
Prop={'w','Kx','Ky','Kxy','p'};     % Properties to add as training data
GP(1)= GP_Plate_Pars(GP(1),Type,AP,Prop);

% Initial value of hyperparameters (D, A, lx, ly, noises)
hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3 10^-3])';

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
Prop={'w','Rx','Kx'};

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

%% Example 2: Curvature-Displacement-based learning - Clamped Boundary Conditions applied
% Reproducibility
rng(10);

% Create boundary condition example
AP.x=[0 AP.x 1];AP.y=[0 AP.y 1];
[AP] = Example2_Analytical_FixedRecPlate(AP);       

% Training parser
Type='Train'; %Training parser

% Boundary conditions
BC.x=[0 AP.b];
BC.y=[0 AP.a];

% Properties to add as training data, including BCs
Prop={'w','Rx','Ry'}; 
GP(2)= GP_Plate_Pars(GP(2),Type,AP,Prop,BC);

% Properties to add as training data, excluding BCs
Prop={'p','Kx','Ky','Kxy'}; 
GP(2)= GP_Plate_Pars(GP(2),Type,AP,Prop); 

% Initial value of hyperparameters (D, A, lx, ly, noises)
hyp=log([1 1 1 exp(1) 10^-3 10^-3 10^-3 10^-3 10^-3 10^-3 10^-3])';

% Define training objective function
fun=@(hyp)GP_Process_Opt(GP(2).Train.X,GP(2).Train.Y,GP(2).Train.f,...
                         GP(2).Kernel,hyp,GP(2).Par.Stabilizer,GP(2).Par.noise);                     

% Standard optimization: gradient descent
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

% Prediction
% Define properties to predict
Type='Pred';
Prop={'w','Rx','Kx'};

% Define locations to predict
Pred.x=linspace(0,AP.a,21); 
Pred.y=linspace(0,AP.b,21);

% Parser - Feed in prediction data
GP(2)= GP_Plate_Pars(GP(2),Type,Pred,Prop);

% Predictions: standard optimizer
[Pred.Mean,Pred.Var] = GP_Process(GP(2).Train.X,GP(2).Train.Y,GP(2).Train.f,...
                                  GP(2).Pred.X,GP(2).Pred.Y,...
                                  GP(2).Kernel,GP(2).Par.hyp,GP(2).Par.Stabilizer,GP(2).Par.noise);
                              
% Predictions: Bayesian optimizer
[Pred.Mean_MCMC,Pred.Var_MCMC] = GP_Process(GP(2).Train.X,GP(2).Train.Y,GP(2).Train.f,...
                                  GP(2).Pred.X,GP(2).Pred.Y,...
                                  GP(2).Kernel,GP(2).Par.hyp_MCMC,GP(2).Par.Stabilizer,GP(2).Par.noise);

% Parser: store results
Type='Results';
GP(2)= GP_Plate_Pars(GP(2),Type,Pred,Prop); 

%% Save results
save('Example2_FixedRecPlate/Example_2.mat', 'GP','AP');

%% Plots
load('Example2_FixedRecPlate/Example_2.mat');
Example2_Plots(GP,AP)
