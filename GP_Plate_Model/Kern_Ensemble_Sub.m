function [Mat] = Kern_Ensemble_Sub(KernFunc,x1,x2,y1,y2,A,L1,L2,D,nu,sigma,BC)
% Evaluates GP kernels on input data

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

Mat=[];
if ~isempty(x1)&&~isempty(x2)
    if nargin(KernFunc)==7
        Mat=KernFunc(A,L1,L2,x1,x2,y1,y2);
    elseif nargin(KernFunc)==8
        Mat=KernFunc(A,D,L1,L2,x1,x2,y1,y2);
    else
        Mat=KernFunc(A,D,L1,L2,nu,x1,x2,y1,y2);
    end
    
    if nargin>10&&sigma~=0%Noise
        NoiseMat=eye(length(x1))*exp(sigma);
        NoiseMat(BC,BC)=0;                  % Boundary conditions
        Mat=Mat+NoiseMat;
    end    
end

end

