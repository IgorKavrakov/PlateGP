function [lik,Der] = GP_Process_Opt(x,y,f,Kernel,hyp,Stabilizer,noise)
% Optimisation function for GP model: get likelihood and possibly kernel
% derivatives given a set of candidate hyperparameters

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

% Copyright (c) Igor Kavrakov, Gledson Rodrigo Tondo, Guido Morgenthal 2025

%This is implementation of GP likelihood by Igor Kavrakov based on
%Rasmussen&Willams Gaussian Processes for Machine Learning (2006) (Algorithm 2.1)

% Reshape the output;
fy=[f.w;f.Rx;f.Ry;f.Kx;f.Ky;f.Kxy;f.p;f.Qx;f.Qy;f.Mx;f.My;f.Mxy];

% Get the kernel value at training points
Kern=Kern_Ensemble(Kernel.Cov,Kernel.BC,Kernel.nu,x,x,y,y,hyp,noise);

% Stabilzer (i.e. jitter term)
Kern=Kern+eye(size(Kern))*Stabilizer;

%Compute cholesky (for inversion - faster). Check Rasmussen
[L_Kern,p]=chol(Kern,'lower');

%Add recovery for matlab-based optimisation if p>0
if p>0; lik=inf; Der=inf; return; end

% Calculate likelihood value
alpha=L_Kern'\(L_Kern\fy);
lik=1/2*fy'*alpha+sum(log(diag((L_Kern))))+length(fy)/2*log(2*pi);

% Check if opmising with derivative-based gradient descent
if nargout==2
    % Base kernel
    Q =  L_Kern'\(L_Kern\eye(size(Kern))) - alpha*alpha';
    Der=hyp*0;
    % A  Derivative
    KernDerA = Kern_Ensemble(Kernel.DerA,Kernel.BC,Kernel.nu,x,x,y,y,hyp,0);       
    Der(1) = trace(Q*KernDerA)/2;
    % LX  Derivative
    KernDerL1= Kern_Ensemble(Kernel.DerL1,Kernel.BC,Kernel.nu,x,x,y,y,hyp,0);      
    Der(2) = trace(Q*KernDerL1)/2;
    % LY  Derivative
    KernDerL2= Kern_Ensemble(Kernel.DerL2,Kernel.BC,Kernel.nu,x,x,y,y,hyp,0);      
    Der(3) = trace(Q*KernDerL2)/2;
    % D  Derivative
    KernDerD = Kern_Ensemble(Kernel.DerD,Kernel.BC,Kernel.nu,x,x,y,y,hyp,0);      
    Der(4) = trace(Q*KernDerD)/2;    
    % Noise derivative
    Der      = Kern_Ensemble_DerNoise(Der,Kernel.BC,x,hyp,Q,noise);   
end
end

