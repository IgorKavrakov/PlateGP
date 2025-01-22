function GP = GP_Plate_Pars(GP,Type,Data,Prop,BC)
% GP Model parser: Organises the GP structure based on the current step:
% pre-training (Init), training, prediction or storing of results.
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

    switch Type
        case 'Init'  % Initialise structures
            % Property list
            Prop={'w','Rx','Ry','Kx','Ky','Kxy','p','Qx','Qy','Mx','My','Mxy'};

            % Assign kernel function
            GP.Kernel=Kern_Exp_Func();

            % Assign Poisson ratio
            GP.Kernel.nu=GP.Par.nu;

            for i=1:length(Prop)
                % Training data: coordinates and quantities
                GP.Train.x.(Prop{i})=[];
                GP.Train.y.(Prop{i})=[];
                GP.Train.X.(Prop{i})=[];
                GP.Train.Y.(Prop{i})=[];
                GP.Train.f.(Prop{i})=[];

                % Prediction data: coordinates and quantities
                GP.Pred.x.(Prop{i})=[];
                GP.Pred.y.(Prop{i})=[];
                GP.Pred.X.(Prop{i})=[];
                GP.Pred.Y.(Prop{i})=[];
                GP.Pred.f.(Prop{i})=[];
                GP.Pred.f_std.(Prop{i})=[];

                % Data re-organised in a meshed format for plotting
                GP.Train.x_mesh.(Prop{i})=[];
                GP.Train.y_mesh.(Prop{i})=[];
                GP.Train.f_mesh.(Prop{i})=[];
                GP.Pred.x_mesh.(Prop{i})=[];
                GP.Pred.y_mesh.(Prop{i})=[];
                GP.Pred.f_mesh.(Prop{i})=[];
                GP.Pred.f_mesh_std.(Prop{i})=[];

                % MCMC parts
                GP.Pred.f_MCMC.(Prop{i})=[];
                GP.Pred.f_std_MCMC.(Prop{i})=[];
                GP.Pred.f_mesh_MCMC.(Prop{i})=[];
                GP.Pred.f_mesh_std_MCMC.(Prop{i})=[];

                % Boundary condition information
                GP.Kernel.BC.(Prop{i})=[];
            end

        case 'Train' % Fill in training data

            % Loop over provided training properties
            for i=1:length(Prop)
                if isfield(Data,Prop{i})&&~isempty(Data.(Prop{i})) 
                    % Training coordinates
                    GP.Train.x.(Prop{i})=Data.x;
                    GP.Train.y.(Prop{i})=Data.y;
                    
                    % Mesh of training coordinates
                    [GP.Train.x_mesh.(Prop{i}),GP.Train.y_mesh.(Prop{i})]=meshgrid(Data.x,Data.y); %Mesh of training coord
                    GP.Train.X.(Prop{i})=GP.Train.x_mesh.(Prop{i})(:); 
                    GP.Train.Y.(Prop{i})=GP.Train.y_mesh.(Prop{i})(:);
                    
                    % Training data
                    GP.Train.f.(Prop{i})=Data.(Prop{i})(:); 

                    % Boundary conditions
                    VBC=zeros(length(GP.Train.x_mesh.(Prop{i})(:)),1); 
                    if nargin>4&&~isempty(BC) 
                        for j=1:length(BC.x) % BCs in x direction
                            VBC(GP.Train.X.(Prop{i})==BC.x(j))=1; 
                        end
                        for j=1:length(BC.y) % BCs in y direction
                            VBC(GP.Train.Y.(Prop{i})==BC.y(j))=1;
                        end
                        GP.Kernel.BC.(Prop{i})=find(VBC==1);  
                    end

                    % Contaminate data additionally with noise (do not
                    % apply noise at BCs)
                    if isfield(Data,'Noise')&&isfield(Data.Noise,Prop{i})
                        Noise=Data.Noise.(Prop{i})(:); 
                        GP.Train.f.(Prop{i})(~VBC)=GP.Train.f.(Prop{i})(~VBC)+Noise(~VBC); 
                    end
                    GP.Train.f_mesh.(Prop{i})=reshape(GP.Train.f.(Prop{i}),size(Data.(Prop{i})));

                elseif  nargin>4&&~isempty(BC)
                    % Mesh of training coordinates
                    [X,Y]=meshgrid(Data.x,Data.y); 
                    X=X(:);Y=Y(:);
                    
                    % Boundary conditions
                    VBC=zeros(length(X),1); 
                    for j=1:length(BC.x) % BCs in x direction
                        VBC(X==BC.x(j))=1; 
                    end
                    for j=1:length(BC.y) % BCs in y direction
                        VBC(Y==BC.y(j))=1; 
                    end
                    
                    % Add boundary conditions
                    GP.Train.Y.(Prop{i})=Y(VBC==1);
                    GP.Train.X.(Prop{i})=X(VBC==1);
                    
                    % Apply zero boundary conditions
                    GP.Train.f.(Prop{i})=GP.Train.X.(Prop{i})*0;                 
                    GP.Kernel.BC.(Prop{i})=(1:length(GP.Train.X.(Prop{i})))';
                end
            end

        case 'Pred' % Fill in prediction data
            % Loop over prediction types
            for i=1:length(Prop)
                % Coordinates for predictions
                GP.Pred.x.(Prop{i})=Data.x(:); 
                GP.Pred.y.(Prop{i})=Data.y(:);

                % Meshing 
                [GP.Pred.x_mesh.(Prop{i}),GP.Pred.y_mesh.(Prop{i})]=meshgrid(GP.Pred.x.(Prop{i}),GP.Pred.y.(Prop{i}));
                GP.Pred.X.(Prop{i})=GP.Pred.x_mesh.(Prop{i})(:); 
                GP.Pred.Y.(Prop{i})=GP.Pred.y_mesh.(Prop{i})(:);
            end

        case 'Results' % Store Results
            % Standard optimisation
            if isfield(Data, 'Mean') && isfield(Data, 'Var')
                Data.STD=sqrt(diag(Data.Var));% Get point standard deviation
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

            % MCMC optimisation
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

