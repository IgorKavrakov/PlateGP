function [Kernel] = Kern_Exp_Func()
% Generates GP kernels based on symbolic differentiation (requires symbolic
% toolbox from MATLAB). By default, kernels are loaded from pre-saved
% Kernel.mat file, for speed.

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


% The kernels & derivatives are already computed using symbolic calculation
% and saved in Kernel.mat
if (3>2)||~license('test','symbolic_toolbox')
    load('Kernel');
else
    % Symbolic computation of the kernels
    % Set up the main kernel
    syms A L1 L2 x1 x2 y1 y2 D nu %A, L1, L2 - hyper parameters. D - Plate parameter. Nu - Poisson. x1,x2,y1,y2 coordinates of point 1 and 2
    Kern=exp(A).*exp(-1/2*(x1-x2)^2./exp(L1)-1/2*(y1-y2)^2./exp(L2)); %Main kernel
    
    % Get separate kernels for the covariance matirx -  Symbolic computation
    % Kernel matrix 12 x 12
    % The superscript D denotes that the kernel includes the flexural rigidity D
    
    % Kern  =[K_ww      K_wRx      K_wRy      K_wKx     K_wKy     K_wKxy    K_wp   K_wQx   K_wQy   K_wMx   K_wMy   K_wMxy;...
    %         K_Rxw     K_RxRx     K_RxRy     K_RxKx    K_RxKy    K_RxKxy   K_Rxp  K_RxQx  K_RxQy  K_RxMx  K_RxMy  K_RxMxy;...
    %         K_Ryw     K_RyRx     K_RyRy     K_RyKx    K_RyKy    K_RyKxy   K_Ryp  K_RyQx  K_RyQy  K_RyMx  K_RyMy  K_RyMxy;...
    %         K_Kxw     K_KxRx     K_KxRy     K_KxKx    K_KxKy    K_KxKxy   K_Kxp    K_KxQx    K_KxQy    K_KxMx    K_KxMy    K_KxMxy;...
    %         K_Kyw     K_KyRx     K_KyRy     K_KyKx    K_KyKy    K_KyKxy   K_Kyp    K_KyQx    K_KyQy    K_KyMx    K_KyMy    K_KyMxy;...
    %         K_Kxyw    K_KxyRx    K_KxyRy    K_KxyKx   K_KxyKy   K_KxyKxy  K_Kxyp   K_KxyQx   K_KxyQy   K_KxyMx   K_KxyMy   K_KxyMxy;...
    %         K_pw      K_pRx      K_pRy      K_pKx     K_pKy     K_pKxy    K_pp     K_pQx     K_pQy     K_pMx     K_pMy     K_pMxy;...
    %         K_Qxw     K_QxRx     K_QxRy     K_QxKx    K_QxKy    K_QxKxy   K_Qxp   K_QxQx    K_QxQy    K_QxMx    K_QxMy    K_QxMxy;...
    %         K_Qyw     K_QyRx     K_QyRy     K_QyKx    K_QyKy    K_QyKxy   K_Qyp    K_QyQx    K_QyQy    K_QyMx    K_QyMy    K_QyMxy;...
    %         K_Mxw     K_MxRx     K_MxRy     K_MxKx    K_MxKy    K_MxKxy   K_Mxp    K_MxQx    K_MxQy    K_MxMx    K_MxMy    K_MxMxy;...
    %         K_Myw     K_MyRx     K_MyRy     K_MyKx    K_MyKy    K_MyKxy   K_Myp    K_MyQx    K_MyQy    K_MyMx    K_MyMy    K_MyMxy;...
    %         K_Mxyw    K_MxyRx    K_MxyRy    K_MxyKx   K_MxyKy   K_MxyKxy  K_Mxyp   K_MxyQx   K_MxyQy   K_MxyMx   K_MxyMy   K_MxyMxy];
    
    
    % Row Displacement w
    Cov.Kern_ww=Kern;
    Cov.Kern_wRx=diff(Kern,x2);
    Cov.Kern_wRy=diff(Kern,y2);
    Cov.Kern_wKx=-diff(Kern,x2,x2);
    Cov.Kern_wKy=-diff(Kern,y2,y2);
    Cov.Kern_wKxy=-2*diff(Kern,x2,y2);
    Cov.Kern_wp=exp(D).*(diff(Kern,x2,x2,x2,x2)+2*diff(Kern,x2,x2,y2,y2)+diff(Kern,y2,y2,y2,y2));
    Cov.Kern_wQx=-exp(D).*(diff(Kern,x2,x2,x2)+diff(Kern,x2,y2,y2));
    Cov.Kern_wQy=-exp(D).*(diff(Kern,x2,x2,y2)+diff(Kern,y2,y2,y2));
    Cov.Kern_wMx=-exp(D).*(diff(Kern,x2,x2)+nu.*diff(Kern,y2,y2));
    Cov.Kern_wMy=-exp(D).*(diff(Kern,y2,y2)+nu.*diff(Kern,x2,x2));
    Cov.Kern_wMxy=exp(D).*(1-nu).*diff(Kern,x2,y2);
    
    % Row Rotation Rx
    Cov.Kern_Rxw=diff(Kern,x1);
    Cov.Kern_RxRx=diff(Kern,x1,x2);
    Cov.Kern_RxRy=diff(Kern,x1,y2);
    Cov.Kern_RxKx=-diff(Kern,x1,x2,x2);
    Cov.Kern_RxKy=-diff(Kern,x1,y2,y2);
    Cov.Kern_RxKxy=-2*diff(Kern,x1,x2,y2);
    Cov.Kern_Rxp=exp(D).*(diff(Kern,x1,x2,x2,x2,x2)+2*diff(Kern,x1,x2,x2,y2,y2)+diff(Kern,x1,y2,y2,y2,y2));
    Cov.Kern_RxQx=-exp(D).*(diff(Kern,x1,x2,x2,x2)+diff(Kern,x1,x2,y2,y2));
    Cov.Kern_RxQy=-exp(D).*(diff(Kern,x1,x2,x2,y2)+diff(Kern,x1,y2,y2,y2));
    Cov.Kern_RxMx=-exp(D).*(diff(Kern,x1,x2,x2)+nu.*diff(Kern,x1,y2,y2));
    Cov.Kern_RxMy=-exp(D).*(diff(Kern,x1,y2,y2)+nu.*diff(Kern,x1,x2,x2));
    Cov.Kern_RxMxy=exp(D).*(1-nu).*diff(Kern,x1,x2,y2);
    
    % Row Rotation Ry
    Cov.Kern_Ryw=diff(Kern,y1);
    Cov.Kern_RyRx=diff(Kern,y1,x2);
    Cov.Kern_RyRy=diff(Kern,y1,y2);
    Cov.Kern_RyKx=-diff(Kern,y1,x2,x2);
    Cov.Kern_RyKy=-diff(Kern,y1,y2,y2);
    Cov.Kern_RyKxy=-2*diff(Kern,y1,x2,y2);
    Cov.Kern_Ryp=exp(D).*(diff(Kern,y1,x2,x2,x2,x2)+2*diff(Kern,y1,x2,x2,y2,y2)+diff(Kern,y1,y2,y2,y2,y2));
    Cov.Kern_RyQx=-exp(D).*(diff(Kern,y1,x2,x2,x2)+diff(Kern,y1,x2,y2,y2));
    Cov.Kern_RyQy=-exp(D).*(diff(Kern,y1,x2,x2,y2)+diff(Kern,y1,y2,y2,y2));
    Cov.Kern_RyMx=-exp(D).*(diff(Kern,y1,x2,x2)+nu.*diff(Kern,y1,y2,y2));
    Cov.Kern_RyMy=-exp(D).*(diff(Kern,y1,y2,y2)+nu.*diff(Kern,y1,x2,x2));
    Cov.Kern_RyMxy=exp(D).*(1-nu).*diff(Kern,y1,x2,y2);
    
    % Row Curvature Kx
    Cov.Kern_Kxw=-diff(Kern,x1,x1);
    Cov.Kern_KxRx=-diff(Kern,x1,x1,x2);
    Cov.Kern_KxRy=-diff(Kern,x1,x1,y2);
    Cov.Kern_KxKx=diff(Kern,x1,x1,x2,x2);
    Cov.Kern_KxKy=diff(Kern,x1,y1,x2,y2);
    Cov.Kern_KxKxy=2*diff(Kern,x1,x1,x2,y2);
    Cov.Kern_Kxp=-exp(D).*(diff(Kern,x1,x1,x2,x2,x2,x2)+2*diff(Kern,x1,x1,x2,x2,y2,y2)+diff(Kern,x1,x1,y2,y2,y2,y2));
    Cov.Kern_KxQx=exp(D).*(diff(Kern,x1,x1,x2,x2,x2)+diff(Kern,x1,x1,x2,y2,y2));
    Cov.Kern_KxQy=exp(D).*(diff(Kern,x1,x1,x2,x2,y2)+diff(Kern,x1,x1,y2,y2,y2));
    Cov.Kern_KxMx=exp(D).*(diff(Kern,x1,x1,x2,x2)+nu.*diff(Kern,x1,x1,y2,y2));
    Cov.Kern_KxMy=exp(D).*(diff(Kern,x1,x1,y2,y2)+nu.*diff(Kern,x1,x1,x2,x2));
    Cov.Kern_KxMxy=-exp(D).*(1-nu).*diff(Kern,x1,x1,x2,y2);
    
    % Row Curvature Ky
    Cov.Kern_Kyw=-diff(Kern,y1,y1);
    Cov.Kern_KyRx=-diff(Kern,y1,y1,x2);
    Cov.Kern_KyRy=-diff(Kern,y1,y1,y2);
    Cov.Kern_KyKx=diff(Kern,x1,y1,x2,y2);
    Cov.Kern_KyKy=diff(Kern,y1,y1,y2,y2);
    Cov.Kern_KyKxy=2*diff(Kern,y1,y1,x2,y2);
    Cov.Kern_Kyp=-exp(D).*(diff(Kern,y1,y1,x2,x2,x2,x2)+2*diff(Kern,y1,y1,x2,x2,y2,y2)+diff(Kern,y1,y1,y2,y2,y2,y2));
    Cov.Kern_KyQx=exp(D).*(diff(Kern,y1,y1,x2,x2,x2)+diff(Kern,y1,y1,x2,y2,y2));
    Cov.Kern_KyQy=exp(D).*(diff(Kern,y1,y1,x2,x2,y2)+diff(Kern,y1,y1,y2,y2,y2));
    Cov.Kern_KyMx=exp(D).*(diff(Kern,y1,y1,x2,x2)+nu.*diff(Kern,y1,y1,y2,y2));
    Cov.Kern_KyMy=exp(D).*(diff(Kern,y1,y1,y2,y2)+nu.*diff(Kern,y1,y1,x2,x2));
    Cov.Kern_KyMxy=-exp(D).*(1-nu).*diff(Kern,y1,y1,x2,y2);
    
    % Row Curvature Kxy
    Cov.Kern_Kxyw=-2*diff(Kern,y1,x1);
    Cov.Kern_KxyRx=-2*diff(Kern,x1,y1,x2);
    Cov.Kern_KxyRy=-2*diff(Kern,x1,y1,y2);
    Cov.Kern_KxyKx=2*diff(Kern,x1,y1,x2,x2);
    Cov.Kern_KxyKy=2*diff(Kern,x1,y1,y2,y2);
    Cov.Kern_KxyKxy=4*diff(Kern,x1,y1,x2,y2);
    Cov.Kern_Kxyp=-2*exp(D).*(diff(Kern,x1,y1,x2,x2,x2,x2)+2*diff(Kern,x1,y1,x2,x2,y2,y2)+diff(Kern,x1,y1,y2,y2,y2,y2));
    Cov.Kern_KxyQx=2*exp(D).*(diff(Kern,x1,y1,x2,x2,x2)+diff(Kern,x1,y1,x2,y2,y2));
    Cov.Kern_KxyQy=2*exp(D).*(diff(Kern,x1,y1,x2,x2,y2)+diff(Kern,x1,y1,y2,y2,y2));
    Cov.Kern_KxyMx=2*exp(D).*(diff(Kern,x1,y1,x2,x2)+nu.*diff(Kern,x1,y1,y2,y2));
    Cov.Kern_KxyMy=2*exp(D).*(diff(Kern,x1,y1,y2,y2)+nu.*diff(Kern,x1,y1,x2,x2));
    Cov.Kern_KxyMxy=-2*exp(D).*(1-nu).*diff(Kern,x1,y1,x2,y2);
    
    % Row load p (in the manuscript - q)
    Cov.Kern_pw=exp(D).*(diff(Kern,x1,x1,x1,x1)+2*diff(Kern,x1,x1,y1,y1)+diff(Kern,y1,y1,y1,y1));
    Cov.Kern_pRx=exp(D).*(diff(Kern,x1,x1,x1,x1,x2)+2*diff(Kern,x1,x1,y1,y1,x2)+diff(Kern,y1,y1,y1,y1,x2));
    Cov.Kern_pRy=exp(D).*(diff(Kern,x1,x1,x1,x1,y2)+2*diff(Kern,x1,x1,y1,y1,y2)+diff(Kern,y1,y1,y1,y1,y2));
    Cov.Kern_pKx=-exp(D).*(diff(Kern,x2,x2,x2,x2,x1,x1)+2*diff(Kern,x1,x1,y1,y1,x2,x2)+diff(Kern,y1,y1,y1,y1,x2,x2));
    Cov.Kern_pKy=-exp(D).*(diff(Kern,x1,x1,x1,x1,y2,y2)+2*diff(Kern,x1,x1,y1,y1,y2,y2)+diff(Kern,y1,y1,y1,y1,y2,y2));
    Cov.Kern_pKxy=-2*exp(D).*(diff(Kern,x1,x1,x1,x1,x2,y2)+2*diff(Kern,x1,x1,y1,y1,x1,y1)+diff(Kern,y1,y1,y1,y1,x2,y2));
    Cov.Kern_pp=(exp(D).^2).*(diff(Kern,x1,x1,x1,x1,x2,x2,x2,x2)+2*diff(Kern,x1,x1,x1,x1,x2,x2,y2,y2)+...
        diff(Kern,x1,x1,x1,x1,y2,y2,y2,y2)+2*diff(Kern,x1,x1,y1,y1,x2,x2,x2,x2)+...
        4*diff(Kern,x1,x1,y1,y1,x2,x2,y2,y2)+2*diff(Kern,x1,x1,y1,y1,y2,y2,y2,y2)+...
        diff(Kern,y1,y1,y1,y1,x2,x2,x2,x2)+2*diff(Kern,y1,y1,y1,y1,x2,x2,y2,y2)+...
        diff(Kern,y1,y1,y1,y1,y2,y2,y2,y2));
    Cov.Kern_pQx=-(exp(D).^2).*(diff(Kern,x1,x1,x1,x1,x2,x2,x2)+diff(Kern,x1,x1,x1,x1,x2,y2,y2)+...
        2*diff(Kern,x1,x1,y1,y1,x2,x2,x2)+2*diff(Kern,x1,x1,y1,y1,x2,y2,y2)+...
        diff(Kern,y1,y1,y1,y1,x2,x2,x2)+diff(Kern,y1,y1,y1,y1,x2,y2,y2));
    Cov.Kern_pQy=-(exp(D).^2).*(diff(Kern,x1,x1,x1,x1,x2,x2,y2)+diff(Kern,x1,x1,x1,x1,y2,y2,y2)+...
        2*diff(Kern,x1,x1,y1,y1,x2,x2,y2)+2*diff(Kern,x1,x1,y1,y1,y2,y2,y2)+...
        diff(Kern,y1,y1,y1,y1,x2,x2,y2)+diff(Kern,y1,y1,y1,y1,y2,y2,y2));
    Cov.Kern_pMx=-(exp(D).^2).*(diff(Kern,x1,x1,x1,x1,x2,x2)+nu.*diff(Kern,x1,x1,x1,x1,y2,y2)+...
        2*diff(Kern,x1,x1,y1,y1,x2,x2)+nu.*2*diff(Kern,x1,x1,y1,y1,y2,y2)+...
        diff(Kern,y1,y1,y1,y1,x2,x2)+nu.*diff(Kern,y1,y1,y1,y1,y2,y2));
    Cov.Kern_pMy=-(exp(D).^2).*(diff(Kern,x1,x1,x1,x1,y2,y2)+nu.*diff(Kern,x1,x1,x1,x1,x2,x2)+...
        2*diff(Kern,x1,x1,y1,y1,y2,y2)+nu.*2*diff(Kern,x1,x1,y1,y1,x2,x2)+...
        diff(Kern,y1,y1,y1,y1,y2,y2)+nu.*diff(Kern,y1,y1,y1,y1,x2,x2));
    Cov.Kern_pMxy=(exp(D).^2).*(1-nu).*(diff(Kern,x1,x1,x1,x1,x2,y2)+2*diff(Kern,x1,x1,y1,y1,x2,y2)+diff(Kern,y1,y1,y1,y1,x2,y2));
    
    % Row shear force Qx
    Cov.Kern_Qxw=-exp(D).*(diff(Kern,x1,x1,x1)+diff(Kern,x1,y1,y1));
    Cov.Kern_QxRx=-exp(D).*(diff(Kern,x1,x1,x1,x2)+diff(Kern,x1,y1,y1,x2));
    Cov.Kern_QxRy=-exp(D).*(diff(Kern,x1,x1,x1,y2)+diff(Kern,x1,y1,y1,y2));
    Cov.Kern_QxKx=exp(D).*(diff(Kern,x1,x1,x1,x2,x2)+diff(Kern,x1,y1,y1,x2,x2));
    Cov.Kern_QxKy=exp(D).*(diff(Kern,x1,x1,x1,y2,y2)+diff(Kern,x1,y1,y1,y2,y2));
    Cov.Kern_QxKxy=2*exp(D).*(diff(Kern,x1,x1,x1,x2,y2)+diff(Kern,x1,y1,y1,x2,y2));
    Cov.Kern_Qxp=-(exp(D).^2).*(diff(Kern,x1,x1,x1,x2,x2,x2,x2)+diff(Kern,x1,y1,y1,x2,x2,x2,x2)+...
        2*diff(Kern,x1,x1,x1,x2,x2,y2,y2)+2*diff(Kern,x1,y1,y1,x2,x2,y2,y2)+...
        diff(Kern,x1,x1,x1,y2,y2,y2,y2)+diff(Kern,x1,y1,y1,y2,y2,y2,y2));
    Cov.Kern_QxQx=(exp(D).^2).*(diff(Kern,x1,x1,x1,x2,x2,x2)+diff(Kern,x1,x1,x1,x2,y2,y2)+...
        diff(Kern,x1,y1,y1,x2,x2,x2)+diff(Kern,x1,y1,y1,x2,y2,y2));
    Cov.Kern_QxQy=(exp(D).^2).*(diff(Kern,x1,x1,x1,x2,x2,y2)+diff(Kern,x1,x1,x1,y2,y2,y2)+...
        diff(Kern,x1,y1,y1,x2,x2,y2)+diff(Kern,x1,y1,y1,y2,y2,y2));
    Cov.Kern_QxMx=(exp(D).^2).*(diff(Kern,x1,x1,x1,x2,x2)+nu.*diff(Kern,x1,x1,x1,y2,y2)+...
        diff(Kern,x1,y1,y1,x2,x2)+nu.*diff(Kern,x1,y1,y1,y2,y2));
    Cov.Kern_QxMy=(exp(D).^2).*(diff(Kern,x1,x1,x1,y2,y2)+nu.*diff(Kern,x1,x1,x1,x2,x2)+...
        diff(Kern,x1,y1,y1,y2,y2)+nu.*diff(Kern,x1,y1,y1,x2,x2));
    Cov.Kern_QxMxy=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,x1,x1,x2,y2)+diff(Kern,x1,y1,y1,x2,y2));
    
    % Row shear force Qy
    Cov.Kern_Qyw=-exp(D).*(diff(Kern,x1,x1,y1)+diff(Kern,y1,y1,y1));
    Cov.Kern_QyRx=-exp(D).*(diff(Kern,x1,x1,y1,x2)+diff(Kern,y1,y1,y1,x2));
    Cov.Kern_QyRy=-exp(D).*(diff(Kern,x1,x1,y1,y2)+diff(Kern,y1,y1,y1,y2));
    Cov.Kern_QyKx=exp(D).*(diff(Kern,x1,x1,y1,x2,x2)+diff(Kern,y1,y1,y1,x2,x2));
    Cov.Kern_QyKy=exp(D).*(diff(Kern,x1,x1,y1,y2,y2)+diff(Kern,y1,y1,y1,y2,y2));
    Cov.Kern_QyKxy=2*exp(D).*(diff(Kern,x1,x1,y1,x2,y2)+diff(Kern,y1,y1,y1,x2,y2));
    Cov.Kern_Qyp=-(exp(D).^2).*(diff(Kern,x1,x1,y1,x2,x2,x2,x2)+diff(Kern,y1,y1,y1,x2,x2,x2,x2)+...
        2*diff(Kern,x1,x1,y1,x2,x2,y2,y2)+2*diff(Kern,y1,y1,y1,x2,x2,y2,y2)+...
        diff(Kern,x1,x1,y1,y2,y2,y2,y2)+diff(Kern,y1,y1,y1,y2,y2,y2,y2));
    Cov.Kern_QyQx=(exp(D).^2).*(diff(Kern,x1,x1,y1,x2,x2,x2)+diff(Kern,y1,y1,y1,x2,x2,x2)+...
        diff(Kern,x1,x1,y1,x2,y2,y2)+diff(Kern,y1,y1,y1,x2,y2,y2));
    Cov.Kern_QyQy=(exp(D).^2).*(diff(Kern,x1,x1,y1,x2,x2,y2)+diff(Kern,x1,x1,y1,y2,y2,y2)+...
        diff(Kern,y1,y1,y1,x2,x2,y2)+diff(Kern,y1,y1,y1,y2,y2,y2));
    Cov.Kern_QyMx=(exp(D).^2).*(diff(Kern,x1,x1,y1,x2,x2)+nu.*diff(Kern,x1,x1,y1,y2,y2)+...
        diff(Kern,y1,y1,y1,x2,x2)+nu.*diff(Kern,y1,y1,y1,y2,y2));
    Cov.Kern_QyMy=(exp(D).^2).*(diff(Kern,x1,x1,y1,y2,y2)+nu.*diff(Kern,x1,x1,y1,x2,x2)+...
        diff(Kern,y1,y1,y1,y2,y2)+nu.*diff(Kern,y1,y1,y1,x2,x2));
    Cov.Kern_QyMxy=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,x1,y1,x2,y2)+diff(Kern,y1,y1,y1,x2,y2));
    % Row Moment Mx
    
    Cov.Kern_Mxw=-exp(D).*(diff(Kern,x1,x1)+nu.*diff(Kern,y1,y1));
    Cov.Kern_MxRx=-exp(D).*(diff(Kern,x1,x1,x2)+nu.*diff(Kern,y1,y1,x2));
    Cov.Kern_MxRy=-exp(D).*(diff(Kern,x1,x1,y2)+nu.*diff(Kern,y1,y1,y2));
    Cov.Kern_MxKx=exp(D).*(diff(Kern,x1,x1,x2,x2)+nu.*diff(Kern,y1,y1,x2,x2));
    Cov.Kern_MxKy=exp(D).*(diff(Kern,x1,x1,y2,y2)+nu.*diff(Kern,y1,y1,y2,y2));
    Cov.Kern_MxKxy=2*exp(D).*(diff(Kern,x1,x1,x2,y2)+nu.*diff(Kern,y1,y1,x2,y2));
    Cov.Kern_Mxp=-(exp(D).^2).*(diff(Kern,x1,x1,x2,x2,x2,x2)+nu.*diff(Kern,y1,y1,x2,x2,x2,x2)+...
        2*diff(Kern,x1,x1,x2,x2,y2,y2)+nu.*2*diff(Kern,y1,y1,x2,x2,y2,y2)+...
        diff(Kern,x1,x1,y2,y2,y2,y2)+nu.*diff(Kern,y1,y1,y2,y2,y2,y2));
    Cov.Kern_MxQx=(exp(D).^2).*(diff(Kern,x1,x1,x2,x2,x2)+nu.*diff(Kern,y1,y1,x2,x2,x2)+...
        diff(Kern,x1,x1,x2,y2,y2)+nu.*diff(Kern,y1,y1,x2,y2,y2));
    Cov.Kern_MxQy=(exp(D).^2).*(diff(Kern,x1,x1,x2,x2,y2)+nu.*diff(Kern,y1,y1,x2,x2,y2)+...
        diff(Kern,x1,x1,y2,y2,y2)+nu.*diff(Kern,y1,y1,y2,y2,y2));
    Cov.Kern_MxMx=(exp(D).^2).*(diff(Kern,x1,x1,x2,x2)+nu.*diff(Kern,x1,x1,y2,y2)+...
        nu*diff(Kern,y1,y1,x2,x2)+nu.^2.*diff(Kern,y1,y1,y2,y2));
    Cov.Kern_MxMy=(exp(D).^2).*(diff(Kern,x1,x1,y2,y2)+nu.*diff(Kern,x1,x1,x2,x2)+...
        nu*diff(Kern,y1,y1,y2,y2)+nu.^2.*diff(Kern,y1,y1,x2,x2));
    Cov.Kern_MxMxy=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,x1,x2,y2)+nu.*diff(Kern,y1,y1,x2,y2));
    
    % Row Moment My
    Cov.Kern_Myw=-exp(D).*(diff(Kern,y1,y1)+nu.*diff(Kern,x1,x1));
    Cov.Kern_MyRx=-exp(D).*(diff(Kern,y1,y1,x2)+nu.*diff(Kern,x1,x1,x2));
    Cov.Kern_MyRy=-exp(D).*(diff(Kern,y1,y1,y2)+nu.*diff(Kern,x1,x1,y2));
    Cov.Kern_MyKx=exp(D).*(diff(Kern,y1,y1,x2,x2)+nu.*diff(Kern,x1,x1,x2,x2));
    Cov.Kern_MyKy=exp(D).*(diff(Kern,y1,y1,y2,y2)+nu.*diff(Kern,x1,x1,y2,y2));
    Cov.Kern_MyKxy=2*exp(D).*(diff(Kern,y1,y1,x2,y2)+nu.*diff(Kern,x1,x1,x2,y2));
    Cov.Kern_Myp=-(exp(D).^2).*(diff(Kern,y1,y1,x2,x2,x2,x2)+nu.*diff(Kern,x1,x1,x2,x2,x2,x2)+...
        2*diff(Kern,y1,y1,x1,x1,y2,y2)+nu.*2*diff(Kern,x1,x1,x2,x2,y2,y2)+...
        diff(Kern,y1,y1,y2,y2,y2,y2)+nu.*diff(Kern,x1,x1,y2,y2,y2,y2));
    Cov.Kern_MyQx=(exp(D).^2).*(diff(Kern,y1,y1,x2,x2,x2)+nu.*diff(Kern,x1,x1,x2,x2,x2)+...
        diff(Kern,y1,y1,x2,y2,y2)+nu.*diff(Kern,x1,x1,x2,y2,y2));
    Cov.Kern_MyQy=(exp(D).^2).*(diff(Kern,y1,y1,x2,x2,y2)+nu.*diff(Kern,x1,x1,x2,x2,y2)+...
        diff(Kern,y1,y1,y2,y2,y2)+nu.*diff(Kern,x1,x1,y2,y2,y2));
    Cov.Kern_MyMx=(exp(D).^2).*(diff(Kern,y1,y1,x2,x2)+nu.*diff(Kern,x1,x1,x2,x2)+...
        nu*diff(Kern,y1,y1,y2,y2)+nu.^2.*diff(Kern,x1,x1,y2,y2));
    Cov.Kern_MyMy=(exp(D).^2).*(diff(Kern,y1,y1,y2,y2)+nu.*diff(Kern,y1,y1,x2,x2)+...
        nu*diff(Kern,x1,x1,y2,y2)+nu.^2.*diff(Kern,x1,x1,x2,x2));
    Cov.Kern_MyMxy=-(exp(D).^2).*(1-nu).*(diff(Kern,y1,y1,x2,y2)+nu.*diff(Kern,x1,x1,x2,y2));
    
    % Row Moment Mxy
    Cov.Kern_Mxyw=exp(D).*(1-nu).*diff(Kern,x1,y1);
    Cov.Kern_MxyRx=exp(D).*(1-nu).*diff(Kern,x1,y1,x1);
    Cov.Kern_MxyRy=exp(D).*(1-nu).*diff(Kern,x1,y1,y2);
    Cov.Kern_MxyKx=-exp(D).*(1-nu).*diff(Kern,x1,y1,x2,x2);
    Cov.Kern_MxyKy=-exp(D).*(1-nu).*diff(Kern,x1,y1,y2,y2);
    Cov.Kern_MxyKxy=-2*exp(D).*(1-nu).*diff(Kern,x1,y1,x2,y2);
    Cov.Kern_Mxyp=(exp(D).^2).*(1-nu).*(diff(Kern,x1,y1,x2,x2,x2,x2)+2*diff(Kern,x1,y1,x2,x2,y2,y2)+diff(Kern,x1,y1,y2,y2,y2,y2));
    Cov.Kern_MxyQx=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,y1,x2,x2,x2)+diff(Kern,x1,y1,x2,y2,y2));
    Cov.Kern_MxyQy=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,y1,x2,x2,y2)+diff(Kern,x1,y1,y2,y2,y2));
    Cov.Kern_MxyMx=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,y1,x2,x2)+nu.*diff(Kern,x1,y1,y2,y2));
    Cov.Kern_MxyMy=-(exp(D).^2).*(1-nu).*(diff(Kern,x1,y1,y2,y2)+nu.*diff(Kern,x1,y1,x2,x2));
    Cov.Kern_MxyMxy=(exp(D).^2).*(1-nu).^2.*diff(Kern,x1,y1,x2,y2);
    
    Kernel.Names=fieldnames(Cov);
    
    % Transfer to function
    for i=1:numel(Kernel.Names)
        Kernel.Cov{i}=matlabFunction(Cov.(Kernel.Names{i}));
    end
    
    % Get separate kernels for the derivative of A for covariance matirx -  Symbolic computation
    for i=1:numel(Kernel.Names)
        DerA{i}=diff(Cov.(Kernel.Names{i}),A);
    end
    
    %Transfer to function
    for i=1:numel(Kernel.Names)
        Kernel.DerA{i}=matlabFunction(DerA{i});
    end
    
    % Get separate kernels for the derivative of L1 for covariance matirx -  Symbolic computation    
    for i=1:numel(Kernel.Names)
        DerL1{i}=diff(Cov.(Kernel.Names{i}),L1);
    end
    
    % Transfer to function
    for i=1:numel(Kernel.Names)
        Kernel.DerL1{i}=matlabFunction(DerL1{i});
    end
    
    % Get separate kernels for the derivative of L2 for covariance matirx -  Symbolic computation
    for i=1:numel(Kernel.Names)
        DerL2{i}=diff(Cov.(Kernel.Names{i}),L2);
    end
    
    % Transfer to function
    for i=1:numel(Kernel.Names)
        Kernel.DerL2{i}=matlabFunction(DerL2{i});
    end
    
    % Get separate kernels for the derivative of D for covariance matirx (Most of them do not contain D) -  Symbolic computation
    KernD={'Kern_wp','Kern_wQx','Kern_wQy','Kern_wMx','Kern_wMy','Kern_wMxy',... %Row 1
        'Kern_Rxp','Kern_RxQx','Kern_RxQy','Kern_RxMx','Kern_RxMy','Kern_RxMxy',... %Row 2
        'Kern_Ryp','Kern_RyQx','Kern_RyQy','Kern_RyMx','Kern_RyMy','Kern_RyMxy',... %Row 3
        'Kern_Kxp','Kern_KxQx','Kern_KxQy','Kern_KxMx','Kern_KxMy','Kern_KxMxy',... %Row 4
        'Kern_Kyp','Kern_KyQx','Kern_KyQy','Kern_KyMx','Kern_KyMy','Kern_KyMxy',... %Row 5
        'Kern_Kxyp','Kern_KxyQx','Kern_KxyQy','Kern_KxyMx','Kern_KxyMy','Kern_KxyMxy',... %Row 6
        'Kern_pw','Kern_pRx','Kern_pRy','Kern_pKx','Kern_pKy','Kern_pKxy','Kern_pp','Kern_pQx','Kern_pQy','Kern_pMx','Kern_pMy','Kern_pMxy',... %Row 7
        'Kern_Qxw','Kern_QxRx','Kern_QxRy','Kern_QxKx','Kern_QxKy','Kern_QxKxy','Kern_Qxyp','Kern_QxQx','Kern_QxQy','Kern_QxMx','Kern_QxMy','Kern_QxMxy',... %Row 8
        'Kern_Qyw','Kern_QyRx','Kern_QyRy','Kern_QyKx','Kern_QyKy','Kern_QyKxy','Kern_Qyp','Kern_QyQx','Kern_QyQy','Kern_QyMx','Kern_QyMy','Kern_QyMxy',... %Row 9
        'Kern_Mxw','Kern_MxRx','Kern_MxRy','Kern_MxKx','Kern_MxKy','Kern_MxKxy','Kern_Mxp','Kern_MxQx','Kern_MxQy','Kern_MxMx','Kern_MxMy','Kern_MxMxy',... %Row 10
        'Kern_Myw','Kern_MyRx','Kern_MyRy','Kern_MyKx','Kern_MyKy','Kern_MyKxy','Kern_Myp','Kern_MyQx','Kern_MyQy','Kern_MyMx','Kern_MyMy','Kern_MyMxy',... %Row 11
        'Kern_Mxyw','Kern_MxyRx','Kern_MxyRy','Kern_MxyKx','Kern_MxyKy','Kern_MxyKxy','Kern_Mxyp','Kern_MxyQx','Kern_MxyQy','Kern_MxyMx','Kern_MxyMy','Kern_MxyMxy'}; %Row 12
    
    DerD=cell(numel(Kernel.Names),1);
    for i=1:numel(Kernel.Names)
        if any(strcmpi(Kernel.Names{i},KernD))
            DerD{i}=diff(Cov.(Kernel.Names{i}),D);
        end
    end
    
    %Transfer to function!
    for i=1:numel(Kernel.Names)
        if any(strcmpi(Kernel.Names{i},KernD))
            Kernel.DerD{i}=matlabFunction(DerD{i});
        else
            Kernel.DerD{i}=@(A,L1,L2,x1,x2,y1,y2)zeros(length(x1),length(x2));
        end
    end
    
    % Save Kernel to file
    Path=regexprep(which('Kern_Exp_Func'),'Kern_Exp_Func.m','');
    save([Path 'Kernel'],'Kernel');
end

end

