function [Der] = Kern_Ensemble_DerNoise(Der,BC,x,hyp,Q,noise)
% Evaluate kernel derivatives w.r.t. the hyperparameters
% Used for gradient based optimisation

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

% Chec for noise
if noise==1
    
    % Count number of inputs
    Nw=length(x.w);
    NRx=length(x.Rx);
    NRy=length(x.Ry);
    NKx=length(x.Kx);
    NKy=length(x.Ky);
    NKxy=length(x.Kxy);
    Np=length(x.p);
    NQx=length(x.Qx);
    NQy=length(x.Qy);
    NMx=length(x.Mx);
    NMy=length(x.My);
    NMxy=length(x.Mxy);
    
    Qdiag=diag(Q);
    
    % Start indice
    Ind=0;n=1;
    
    % Displacements
    if Nw&&(Nw~=length(BC.w))
        Qder=Qdiag(1+Ind:Ind+Nw);
        Qder(BC.w)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+Nw;
    
    % X rotations
    if NRx&&(NRx~=length(BC.Rx))
        Qder=Qdiag(1+Ind:Ind+NRx);
        Qder(BC.Rx)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+NRx;
    
    % Y Rotations
    if NRy&&(NRy~=length(BC.Ry))
        Qder=Qdiag(1+Ind:Ind+NRy);
        Qder(BC.Ry)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+NRy;
    
    % X curvature
    if NKx&&(NKx~=length(BC.Kx))
        Qder=Qdiag(1+Ind:Ind+NKx);
        Qder(BC.Kx)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+NKx;
    
    % Y curvature
    if NKy&&(NKy~=length(BC.Ky))
        Qder=Qdiag(1+Ind:Ind+NKy);
        Qder(BC.Ky)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+NKy;
    
    % XY curvature
    if NKxy&&(NKxy~=length(BC.Kxy))
        Qder=Qdiag(1+Ind:Ind+NKxy);
        Qder(BC.Kxy)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
        n=n+1;
    end
    
    Ind=Ind+NKxy;
    
    % Load
    if Np&&(Np~=length(BC.p))
        Qder=Qdiag(1+Ind:Ind+Np);
        Qder(BC.p)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
    
    Ind=Ind+Np;
    
    % X shear
    if NQx&&(NQx~=length(BC.Qx))
        Qder=Qdiag(1+Ind:Ind+NQx);
        Qder(BC.Qx)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
    
    Ind=Ind+NQx;
    
    % Y shear
    if NQy&&(NQy~=length(BC.Qy))
        Qder=Qdiag(1+Ind:Ind+NQy);
        Qder(BC.Qy)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
    
    Ind=Ind+NQy;
    
    % X moment
    if NMx&&(Nx~=length(BC.Mx))
        Qder=Qdiag(1+Ind:Ind+Nx);
        Qder(BC.Mx)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
    
    Ind=Ind+NMx;
    
    % Y moment
    if NMy&&(NMy~=length(BC.My))
        Qder=Qdiag(1+Ind:Ind+NMy);
        Qder(BC.My)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
    
    Ind=Ind+NMy;
    
    % XY moment
    if NMxy&&(NMxy~=length(BC.Mxy))
        Qder=Qdiag(1+Ind:Ind+NMxyMxy);
        Qder(BC.Mxy)=[];
        Der(4+n)=exp(hyp(4+n))*sum(Qder)/2;
    end
end

end


