function GP = GP_Plate_Pars(GP,Type,Data,Prop,BC)
% GP Model parser

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

switch Type
    case 'Init'  % Initialise structures
    Prop={'w','Rx','Ry','Kx','Ky','Kxy','p','Qx','Qy','Mx','My','Mxy'};
    GP.Kernel=Kern_Exp_Func();
    GP.Kernel.nu=GP.Par.nu; % Place the known Poisson coeff. in the Kernel
    
    for i=1:length(Prop)
       GP.Train.x.(Prop{i})=[]; %Unique coordinate for training
       GP.Train.y.(Prop{i})=[];       
       GP.Train.X.(Prop{i})=[]; %List of coordinates fot training
       GP.Train.Y.(Prop{i})=[];
       GP.Train.f.(Prop{i})=[]; %List of prediction values
       
       GP.Pred.x.(Prop{i})=[]; %Unique coordinate for training
       GP.Pred.y.(Prop{i})=[];       
       GP.Pred.X.(Prop{i})=[]; %List of coordinates fot training
       GP.Pred.Y.(Prop{i})=[];
       GP.Pred.f.(Prop{i})=[]; %List of prediction values         
       GP.Pred.f_std.(Prop{i})=[]; %Standard deviation

       %Plotting data org (Redundant, but conveniant)
       GP.Train.x_mesh.(Prop{i})=[]; 
       GP.Train.y_mesh.(Prop{i})=[];
       GP.Train.f_mesh.(Prop{i})=[]; 
       
       GP.Pred.x_mesh.(Prop{i})=[]; %Mesh points for plotting -> Data reorganised
       GP.Pred.y_mesh.(Prop{i})=[];
       GP.Pred.f_mesh.(Prop{i})=[];
       GP.Pred.f_mesh_std.(Prop{i})=[];
       
       % MCMC parts
       GP.Pred.f_MCMC.(Prop{i})=[];
       GP.Pred.f_std_MCMC.(Prop{i})=[];
       GP.Pred.f_mesh_MCMC.(Prop{i})=[];
       GP.Pred.f_mesh_std_MCMC.(Prop{i})=[];
       
       GP.Kernel.BC.(Prop{i})=[]; %BC to be assigned to the kernel (i.e. GP model)
    end
    
    case 'Train' %% Fill in training data

     for i=1:length(Prop)
       if isfield(Data,Prop{i})&&~isempty(Data.(Prop{i})) % There is data for the selected Prop 
           GP.Train.x.(Prop{i})=Data.x; %Unique coordinate for training
           GP.Train.y.(Prop{i})=Data.y;

           [GP.Train.x_mesh.(Prop{i}),GP.Train.y_mesh.(Prop{i})]=meshgrid(Data.x,Data.y); %Mesh of training coord
           GP.Train.X.(Prop{i})=GP.Train.x_mesh.(Prop{i})(:); %List of coordinates fot training
           GP.Train.Y.(Prop{i})=GP.Train.y_mesh.(Prop{i})(:);
           GP.Train.f.(Prop{i})=Data.(Prop{i})(:); %List of prediction values   

           %Boundary conditions
           VBC=zeros(length(GP.Train.x_mesh.(Prop{i})(:)),1); %Vector containing BC 
           if nargin>4&&~isempty(BC) %Find boundary conditions indices if they are active
               for j=1:length(BC.x)
                 VBC(GP.Train.X.(Prop{i})==BC.x(j))=1; % Boundary conditions in x direction along BC.x(j,1)
               end
               for j=1:length(BC.y)
                 VBC(GP.Train.Y.(Prop{i})==BC.y(j))=1;
               end       
               GP.Kernel.BC.(Prop{i})=find(VBC==1);  %Find indices of active boundary conditions
           end

           %Contaminate data additionally with noise? - if there are no BC!
           if isfield(Data,'Noise')&&isfield(Data.Noise,Prop{i}) % Add noise if needed (each variable has specific noise
               Noise=Data.Noise.(Prop{i})(:); %Separate noise as read
               GP.Train.f.(Prop{i})(~VBC)=GP.Train.f.(Prop{i})(~VBC)+Noise(~VBC); %Apply only where no BC
           end       
           GP.Train.f_mesh.(Prop{i})=reshape(GP.Train.f.(Prop{i}),size(Data.(Prop{i})));
           
       elseif  nargin>4&&~isempty(BC) %Apply zero boundary conditions for the selected property    
           [X,Y]=meshgrid(Data.x,Data.y); %Mesh of training coord 
           X=X(:);Y=Y(:);
           
           VBC=zeros(length(X),1); %Vector containing BC            
           for j=1:length(BC.x)
              VBC(X==BC.x(j))=1; % Boundary conditions in x direction along BC.x(j,1)
           end 
           for j=1:length(BC.y)
              VBC(Y==BC.y(j))=1; % Boundary conditions in x direction along BC.x(j,1)
           end 
           GP.Train.Y.(Prop{i})=Y(VBC==1);
           GP.Train.X.(Prop{i})=X(VBC==1);        
           
           GP.Train.f.(Prop{i})=GP.Train.X.(Prop{i})*0;                 %Apply zero boundary conditions               
           GP.Kernel.BC.(Prop{i})=(1:length(GP.Train.X.(Prop{i})))';  %Find indices of active boundary conditions
       end      
     end
     
    case 'Pred' %% Fill in prediction data
     
     for i=1:length(Prop)
       GP.Pred.x.(Prop{i})=Data.x(:); %Unique coordinate for training
       GP.Pred.y.(Prop{i})=Data.y(:);

       [GP.Pred.x_mesh.(Prop{i}),GP.Pred.y_mesh.(Prop{i})]=meshgrid(GP.Pred.x.(Prop{i}),GP.Pred.y.(Prop{i}));
       GP.Pred.X.(Prop{i})=GP.Pred.x_mesh.(Prop{i})(:); %List of coordinates fot training
       GP.Pred.Y.(Prop{i})=GP.Pred.y_mesh.(Prop{i})(:);
     end
      
    case 'Results' % Store Results
     if isfield(Data, 'Mean') && isfield(Data, 'Var')
         Data.STD=sqrt(diag(Data.Var));%Get point standard deviation
         N=0;
         for i=1:length(Prop)
           Nx=length(GP.Pred.x.(Prop{i}));
           Ny=length(GP.Pred.y.(Prop{i}));
           GP.Pred.f.(Prop{i})=Data.Mean(1+N:N+Nx*Ny); %Store mean prediction
           GP.Pred.f_std.(Prop{i})=Data.STD(1+N:N+Nx*Ny); %Store standard deviation     
           GP.Pred.f_mesh.(Prop{i})=reshape(GP.Pred.f.(Prop{i}),[Nx,Ny]); % Store mean prediction - mesh for plotting
           GP.Pred.f_mesh_std.(Prop{i})=reshape(GP.Pred.f_std.(Prop{i}),[Nx,Ny]); % Store mean prediction - mesh for plotting

           N=N+Nx*Ny;
         end      
     end
     
     if isfield(Data, 'Mean_MCMC') && isfield(Data, 'Var_MCMC')
         Data.STD_MCMC=sqrt(diag(Data.Var_MCMC));%Get point standard deviation
         N=0;
         for i=1:length(Prop)
           Nx=length(GP.Pred.x.(Prop{i}));
           Ny=length(GP.Pred.y.(Prop{i}));
           GP.Pred.f_MCMC.(Prop{i})=Data.Mean_MCMC(1+N:N+Nx*Ny); %Store mean prediction
           GP.Pred.f_std_MCMC.(Prop{i})=Data.STD_MCMC(1+N:N+Nx*Ny); %Store standard deviation     
           GP.Pred.f_mesh_MCMC.(Prop{i})=reshape(GP.Pred.f_MCMC.(Prop{i}),[Nx,Ny]); % Store mean prediction - mesh for plotting
           GP.Pred.f_mesh_std_MCMC.(Prop{i})=reshape(GP.Pred.f_std_MCMC.(Prop{i}),[Nx,Ny]); % Store mean prediction - mesh for plotting

           N=N+Nx*Ny;
         end    
     end
        
end
end

