function [Kern] = Kern_Ensemble(Kernel,BC,nu,x1,x2,y1,y2,hyp,noise)
% Wrapper function for kernel evaluation

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

%hyp=exp(hyp);
A=hyp(1);L1=hyp(2);L2=hyp(3);D=hyp(4);
if noise==1 %Identify noise?
   n=1;
   if ~isempty(x1.w)&&(length(BC.w)~=length(x1.w));      sigma_w=hyp(4+n);  n=n+1; else sigma_w=0;   end
   if ~isempty(x1.Rx)&&(length(BC.Rx)~=length(x1.Rx));   sigma_Rx=hyp(4+n); n=n+1; else sigma_Rx=0;  end
   if ~isempty(x1.Ry)&&(length(BC.Ry)~=length(x1.Ry));   sigma_Ry=hyp(4+n); n=n+1; else sigma_Ry=0;  end   
   if ~isempty(x1.Kx)&&(length(BC.Kx)~=length(x1.Kx));   sigma_Kx=hyp(4+n); n=n+1; else sigma_Kx=0;  end
   if ~isempty(x1.Ky)&&(length(BC.Ky)~=length(x1.Ky));   sigma_Ky=hyp(4+n); n=n+1; else sigma_Ky=0;  end
   if ~isempty(x1.Kxy)&&(length(BC.Kxy)~=length(x1.Kxy));sigma_Kxy=hyp(4+n);n=n+1; else sigma_Kxy=0; end
   if ~isempty(x1.p)&&(length(BC.p)~=length(x1.p));      sigma_p=hyp(4+n);  n=n+1; else sigma_p=0;   end
   if ~isempty(x1.Qx)&&(length(BC.Qx)~=length(x1.p));    sigma_Qx=hyp(4+n); n=n+1; else sigma_Qx=0;  end
   if ~isempty(x1.Qy)&&(length(BC.Qy)~=length(x1.p));    sigma_Qy=hyp(4+n); n=n+1; else sigma_Qy=0;  end
   if ~isempty(x1.Mx)&&(length(BC.Mx)~=length(x1.p));    sigma_Mx=hyp(4+n); n=n+1; else sigma_Mx=0;  end
   if ~isempty(x1.My)&&(length(BC.My)~=length(x1.p));    sigma_My=hyp(4+n); n=n+1; else sigma_My=0;  end
   if ~isempty(x1.Mxy)&&(length(BC.Mxy)~=length(x1.p));  sigma_Mxy=hyp(4+n);       else sigma_Mxy=0; end
else
   sigma_w=0;sigma_Rx=0;sigma_Ry=0;sigma_Kx=0;sigma_Ky=0;sigma_Kxy=0;sigma_p=0;sigma_Qx=0;sigma_Qy=0;sigma_Mx=0;sigma_My=0;sigma_Mxy=0;
end

if isempty(BC)
   BC.w=[]; BC.Rx=[]; BC.Ry=[]; BC.Kx=[]; BC.Ky=[]; BC.Kxy=[]; BC.p=[]; BC.Mx=[]; BC.My=[]; BC.Mxy=[]; BC.Qx=[]; BC.Qy=[];
end
%Get contributions
% Row Displacement w
K_ww    =Kern_Ensemble_Sub(Kernel{1}  ,x1.w,  x2.w',  y1.w,  y2.w',  A,L1,L2,[],[],sigma_w,BC.w);
K_wRx   =Kern_Ensemble_Sub(Kernel{2}  ,x1.w,  x2.Rx', y1.w,  y2.Rx', A,L1,L2);
K_wRy   =Kern_Ensemble_Sub(Kernel{3}  ,x1.w,  x2.Ry', y1.w,  y2.Ry', A,L1,L2);
K_wKx   =Kern_Ensemble_Sub(Kernel{4}  ,x1.w,  x2.Kx', y1.w,  y2.Kx', A,L1,L2);
K_wKy   =Kern_Ensemble_Sub(Kernel{5}  ,x1.w,  x2.Ky', y1.w,  y2.Ky', A,L1,L2);
K_wKxy  =Kern_Ensemble_Sub(Kernel{6}  ,x1.w,  x2.Kxy',y1.w,  y2.Kxy',A,L1,L2);
K_wp    =Kern_Ensemble_Sub(Kernel{7}  ,x1.w,  x2.p'  ,y1.w,  y2.p',  A,L1,L2,D);
K_wQx   =Kern_Ensemble_Sub(Kernel{8}  ,x1.w,  x2.Qx' ,y1.w,  y2.Qx', A,L1,L2,D);
K_wQy   =Kern_Ensemble_Sub(Kernel{9}  ,x1.w,  x2.Qy' ,y1.w,  y2.Qy', A,L1,L2,D);
K_wMx   =Kern_Ensemble_Sub(Kernel{10} ,x1.w,  x2.Mx' ,y1.w,  y2.Mx', A,L1,L2,D,nu);
K_wMy   =Kern_Ensemble_Sub(Kernel{11} ,x1.w,  x2.My' ,y1.w,  y2.My', A,L1,L2,D,nu);
K_wMxy  =Kern_Ensemble_Sub(Kernel{12} ,x1.w,  x2.Mxy',y1.w,  y2.Mxy',A,L1,L2,D,nu);

% Row Rotation Rx 
K_Rxw   =Kern_Ensemble_Sub(Kernel{13} ,x1.Rx, x2.w',  y1.Rx, y2.w',  A,L1,L2);
K_RxRx  =Kern_Ensemble_Sub(Kernel{14} ,x1.Rx, x2.Rx', y1.Rx, y2.Rx', A,L1,L2,[],[],sigma_Rx,BC.Rx);
K_RxRy  =Kern_Ensemble_Sub(Kernel{15} ,x1.Rx, x2.Ry', y1.Rx, y2.Ry', A,L1,L2);
K_RxKx  =Kern_Ensemble_Sub(Kernel{16} ,x1.Rx, x2.Kx', y1.Rx, y2.Kx', A,L1,L2);
K_RxKy  =Kern_Ensemble_Sub(Kernel{17} ,x1.Rx, x2.Ky', y1.Rx, y2.Ky', A,L1,L2);
K_RxKxy =Kern_Ensemble_Sub(Kernel{18} ,x1.Rx, x2.Kxy',y1.Rx, y2.Kxy',A,L1,L2);
K_Rxp   =Kern_Ensemble_Sub(Kernel{19} ,x1.Rx, x2.p'  ,y1.Rx, y2.p',  A,L1,L2,D);
K_RxQx  =Kern_Ensemble_Sub(Kernel{20} ,x1.Rx, x2.Qx' ,y1.Rx, y2.Qx', A,L1,L2,D);
K_RxQy  =Kern_Ensemble_Sub(Kernel{21} ,x1.Rx, x2.Qy' ,y1.Rx, y2.Qy', A,L1,L2,D);
K_RxMx  =Kern_Ensemble_Sub(Kernel{22} ,x1.Rx, x2.Mx' ,y1.Rx, y2.Mx', A,L1,L2,D,nu);
K_RxMy  =Kern_Ensemble_Sub(Kernel{23} ,x1.Rx, x2.My' ,y1.Rx, y2.My', A,L1,L2,D,nu);
K_RxMxy =Kern_Ensemble_Sub(Kernel{24} ,x1.Rx, x2.Mxy',y1.Rx, y2.Mxy',A,L1,L2,D,nu);

% Row Rotation Ry 
K_Ryw   =Kern_Ensemble_Sub(Kernel{25} ,x1.Ry, x2.w',  y1.Ry, y2.w',  A,L1,L2);
K_RyRx  =Kern_Ensemble_Sub(Kernel{26} ,x1.Ry, x2.Rx', y1.Ry, y2.Rx', A,L1,L2);
K_RyRy  =Kern_Ensemble_Sub(Kernel{27} ,x1.Ry, x2.Ry', y1.Ry, y2.Ry', A,L1,L2,[],[],sigma_Ry,BC.Ry);
K_RyKx  =Kern_Ensemble_Sub(Kernel{28} ,x1.Ry, x2.Kx', y1.Ry, y2.Kx', A,L1,L2);
K_RyKy  =Kern_Ensemble_Sub(Kernel{29} ,x1.Ry, x2.Ky', y1.Ry, y2.Ky', A,L1,L2);
K_RyKxy =Kern_Ensemble_Sub(Kernel{30} ,x1.Ry, x2.Kxy',y1.Ry, y2.Kxy',A,L1,L2);
K_Ryp   =Kern_Ensemble_Sub(Kernel{31} ,x1.Ry, x2.p'  ,y1.Ry, y2.p',  A,L1,L2,D);
K_RyQx  =Kern_Ensemble_Sub(Kernel{32} ,x1.Ry, x2.Qx' ,y1.Ry, y2.Qx', A,L1,L2,D);
K_RyQy  =Kern_Ensemble_Sub(Kernel{33} ,x1.Ry, x2.Qy' ,y1.Ry, y2.Qy', A,L1,L2,D);
K_RyMx  =Kern_Ensemble_Sub(Kernel{34} ,x1.Ry, x2.Mx' ,y1.Ry, y2.Mx', A,L1,L2,D,nu);
K_RyMy  =Kern_Ensemble_Sub(Kernel{35} ,x1.Ry, x2.My' ,y1.Ry, y2.My', A,L1,L2,D,nu);
K_RyMxy =Kern_Ensemble_Sub(Kernel{36} ,x1.Ry, x2.Mxy',y1.Ry, y2.Mxy',A,L1,L2,D,nu);

% Row Curvature Kx 
K_Kxw   =Kern_Ensemble_Sub(Kernel{37} ,x1.Kx, x2.w',  y1.Kx, y2.w',  A,L1,L2);
K_KxRx  =Kern_Ensemble_Sub(Kernel{38} ,x1.Kx, x2.Rx', y1.Kx, y2.Rx', A,L1,L2);
K_KxRy  =Kern_Ensemble_Sub(Kernel{39} ,x1.Kx, x2.Ry', y1.Kx, y2.Ry', A,L1,L2);
K_KxKx  =Kern_Ensemble_Sub(Kernel{40} ,x1.Kx, x2.Kx', y1.Kx, y2.Kx', A,L1,L2,[],[],sigma_Kx,BC.Kx);
K_KxKy  =Kern_Ensemble_Sub(Kernel{41} ,x1.Kx, x2.Ky', y1.Kx, y2.Ky', A,L1,L2);
K_KxKxy =Kern_Ensemble_Sub(Kernel{42} ,x1.Kx, x2.Kxy',y1.Kx, y2.Kxy',A,L1,L2);
K_Kxp   =Kern_Ensemble_Sub(Kernel{43} ,x1.Kx, x2.p'  ,y1.Kx, y2.p',  A,L1,L2,D);
K_KxQx  =Kern_Ensemble_Sub(Kernel{44} ,x1.Kx, x2.Qx' ,y1.Kx, y2.Qx', A,L1,L2,D);
K_KxQy  =Kern_Ensemble_Sub(Kernel{45} ,x1.Kx, x2.Qy' ,y1.Kx, y2.Qy', A,L1,L2,D);
K_KxMx  =Kern_Ensemble_Sub(Kernel{46} ,x1.Kx, x2.Mx' ,y1.Kx, y2.Mx', A,L1,L2,D,nu);
K_KxMy  =Kern_Ensemble_Sub(Kernel{47} ,x1.Kx, x2.My' ,y1.Kx, y2.My', A,L1,L2,D,nu);
K_KxMxy =Kern_Ensemble_Sub(Kernel{48} ,x1.Kx, x2.Mxy',y1.Kx, y2.Mxy',A,L1,L2,D,nu);

% Row Curvature Ky 
K_Kyw   =Kern_Ensemble_Sub(Kernel{49} ,x1.Ky, x2.w',  y1.Ky, y2.w',  A,L1,L2);
K_KyRx  =Kern_Ensemble_Sub(Kernel{50} ,x1.Ky, x2.Rx', y1.Ky, y2.Rx', A,L1,L2);
K_KyRy  =Kern_Ensemble_Sub(Kernel{51} ,x1.Ky, x2.Ry', y1.Ky, y2.Ry', A,L1,L2);
K_KyKx  =Kern_Ensemble_Sub(Kernel{52} ,x1.Ky, x2.Kx', y1.Ky, y2.Kx', A,L1,L2);
K_KyKy  =Kern_Ensemble_Sub(Kernel{53} ,x1.Ky, x2.Ky', y1.Ky, y2.Ky', A,L1,L2,[],[],sigma_Ky,BC.Ky);
K_KyKxy =Kern_Ensemble_Sub(Kernel{54} ,x1.Ky, x2.Kxy',y1.Ky, y2.Kxy',A,L1,L2);
K_Kyp   =Kern_Ensemble_Sub(Kernel{55} ,x1.Ky, x2.p'  ,y1.Ky, y2.p',  A,L1,L2,D);
K_KyQx  =Kern_Ensemble_Sub(Kernel{56} ,x1.Ky, x2.Qx' ,y1.Ky, y2.Qx', A,L1,L2,D);
K_KyQy  =Kern_Ensemble_Sub(Kernel{57} ,x1.Ky, x2.Qy' ,y1.Ky, y2.Qy', A,L1,L2,D);
K_KyMx  =Kern_Ensemble_Sub(Kernel{58} ,x1.Ky, x2.Mx' ,y1.Ky, y2.Mx', A,L1,L2,D,nu);
K_KyMy  =Kern_Ensemble_Sub(Kernel{59} ,x1.Ky, x2.My' ,y1.Ky, y2.My', A,L1,L2,D,nu);
K_KyMxy =Kern_Ensemble_Sub(Kernel{60} ,x1.Ky, x2.Mxy',y1.Ky, y2.Mxy',A,L1,L2,D,nu);

% Row Curvature Kxy 
K_Kxyw  =Kern_Ensemble_Sub(Kernel{61} ,x1.Kxy,x2.w',  y1.Kxy,y2.w',  A,L1,L2);
K_KxyRx =Kern_Ensemble_Sub(Kernel{62} ,x1.Kxy,x2.Rx', y1.Kxy,y2.Rx', A,L1,L2);
K_KxyRy =Kern_Ensemble_Sub(Kernel{63} ,x1.Kxy,x2.Ry', y1.Kxy,y2.Ry', A,L1,L2);
K_KxyKx =Kern_Ensemble_Sub(Kernel{64} ,x1.Kxy,x2.Kx', y1.Kxy,y2.Kx', A,L1,L2);
K_KxyKy =Kern_Ensemble_Sub(Kernel{65} ,x1.Kxy,x2.Ky', y1.Kxy,y2.Ky', A,L1,L2);
K_KxyKxy=Kern_Ensemble_Sub(Kernel{66} ,x1.Kxy,x2.Kxy',y1.Kxy,y2.Kxy',A,L1,L2,[],[],sigma_Kxy,BC.Kxy);
K_Kxyp  =Kern_Ensemble_Sub(Kernel{67} ,x1.Kxy,x2.p'  ,y1.Kxy,y2.p',  A,L1,L2,D);
K_KxyQx =Kern_Ensemble_Sub(Kernel{68} ,x1.Kxy,x2.Qx' ,y1.Kxy,y2.Qx', A,L1,L2,D);
K_KxyQy =Kern_Ensemble_Sub(Kernel{69} ,x1.Kxy,x2.Qy' ,y1.Kxy,y2.Qy', A,L1,L2,D);
K_KxyMx =Kern_Ensemble_Sub(Kernel{70} ,x1.Kxy,x2.Mx' ,y1.Kxy,y2.Mx', A,L1,L2,D,nu);
K_KxyMy =Kern_Ensemble_Sub(Kernel{71} ,x1.Kxy,x2.My' ,y1.Kxy,y2.My', A,L1,L2,D,nu);
K_KxyMxy=Kern_Ensemble_Sub(Kernel{72} ,x1.Kxy,x2.Mxy',y1.Kxy,y2.Mxy',A,L1,L2,D,nu);

% Row load p (in the manuscript - q)
K_pw    =Kern_Ensemble_Sub(Kernel{73} ,x1.p,  x2.w',  y1.p,  y2.w',  A,L1,L2,D);
K_pRx   =Kern_Ensemble_Sub(Kernel{74} ,x1.p,  x2.Rx', y1.p,  y2.Rx', A,L1,L2,D);
K_pRy   =Kern_Ensemble_Sub(Kernel{75} ,x1.p,  x2.Ry', y1.p,  y2.Ry', A,L1,L2,D);
K_pKx   =Kern_Ensemble_Sub(Kernel{76} ,x1.p,  x2.Kx', y1.p,  y2.Kx', A,L1,L2,D);
K_pKy   =Kern_Ensemble_Sub(Kernel{77} ,x1.p,  x2.Ky', y1.p,  y2.Ky', A,L1,L2,D);
K_pKxy  =Kern_Ensemble_Sub(Kernel{78} ,x1.p,  x2.Kxy',y1.p,  y2.Kxy',A,L1,L2,D);
K_pp    =Kern_Ensemble_Sub(Kernel{79} ,x1.p,  x2.p'  ,y1.p,  y2.p',  A,L1,L2,D,[],sigma_p,[]);
K_pQx   =Kern_Ensemble_Sub(Kernel{80} ,x1.p,  x2.Qx' ,y1.p,  y2.Qx', A,L1,L2,D);
K_pQy   =Kern_Ensemble_Sub(Kernel{81} ,x1.p,  x2.Qy' ,y1.p,  y2.Qy', A,L1,L2,D);
K_pMx   =Kern_Ensemble_Sub(Kernel{82} ,x1.p,  x2.Mx' ,y1.p,  y2.Mx', A,L1,L2,D,nu);
K_pMy   =Kern_Ensemble_Sub(Kernel{83} ,x1.p,  x2.My' ,y1.p,  y2.My', A,L1,L2,D,nu);
K_pMxy  =Kern_Ensemble_Sub(Kernel{84} ,x1.p,  x2.Mxy',y1.p,  y2.Mxy',A,L1,L2,D,nu);

% Row shear force Qx 
K_Qxw   =Kern_Ensemble_Sub(Kernel{85} ,x1.Qx, x2.w',  y1.Qx, y2.w',  A,L1,L2,D);
K_QxRx  =Kern_Ensemble_Sub(Kernel{86} ,x1.Qx, x2.Rx', y1.Qx, y2.Rx', A,L1,L2,D);
K_QxRy  =Kern_Ensemble_Sub(Kernel{87} ,x1.Qx, x2.Ry', y1.Qx, y2.Ry', A,L1,L2,D);
K_QxKx  =Kern_Ensemble_Sub(Kernel{88} ,x1.Qx, x2.Kx', y1.Qx, y2.Kx', A,L1,L2,D);
K_QxKy  =Kern_Ensemble_Sub(Kernel{89} ,x1.Qx, x2.Ky', y1.Qx, y2.Ky', A,L1,L2,D);
K_QxKxy =Kern_Ensemble_Sub(Kernel{90} ,x1.Qx, x2.Kxy',y1.Qx, y2.Kxy',A,L1,L2,D);
K_Qxp   =Kern_Ensemble_Sub(Kernel{91} ,x1.Qx, x2.p'  ,y1.Qx, y2.p',  A,L1,L2,D);
K_QxQx  =Kern_Ensemble_Sub(Kernel{92} ,x1.Qx, x2.Qx' ,y1.Qx, y2.Qx', A,L1,L2,D,[],sigma_Qx,BC.Qx);
K_QxQy  =Kern_Ensemble_Sub(Kernel{93} ,x1.Qx, x2.Qy' ,y1.Qx, y2.Qy', A,L1,L2,D);
K_QxMx  =Kern_Ensemble_Sub(Kernel{94} ,x1.Qx, x2.Mx' ,y1.Qx, y2.Mx', A,L1,L2,D,nu);
K_QxMy  =Kern_Ensemble_Sub(Kernel{95} ,x1.Qx, x2.My' ,y1.Qx, y2.My', A,L1,L2,D,nu);
K_QxMxy =Kern_Ensemble_Sub(Kernel{96} ,x1.Qx, x2.Mxy',y1.Qx, y2.Mxy',A,L1,L2,D,nu);

% Row shear force Qy 
K_Qyw   =Kern_Ensemble_Sub(Kernel{97} ,x1.Qy, x2.w',  y1.Qy, y2.w',  A,L1,L2,D);
K_QyRx  =Kern_Ensemble_Sub(Kernel{98} ,x1.Qy, x2.Rx', y1.Qy, y2.Rx', A,L1,L2,D);
K_QyRy  =Kern_Ensemble_Sub(Kernel{99} ,x1.Qy, x2.Ry', y1.Qy, y2.Ry', A,L1,L2,D);
K_QyKx  =Kern_Ensemble_Sub(Kernel{100},x1.Qy, x2.Kx', y1.Qy, y2.Kx', A,L1,L2,D);
K_QyKy  =Kern_Ensemble_Sub(Kernel{101},x1.Qy, x2.Ky', y1.Qy, y2.Ky', A,L1,L2,D);
K_QyKxy =Kern_Ensemble_Sub(Kernel{102},x1.Qy, x2.Kxy',y1.Qy, y2.Kxy',A,L1,L2,D);
K_Qyp   =Kern_Ensemble_Sub(Kernel{103},x1.Qy, x2.p'  ,y1.Qy, y2.p',  A,L1,L2,D);
K_QyQx  =Kern_Ensemble_Sub(Kernel{104},x1.Qy, x2.Qx' ,y1.Qy, y2.Qx', A,L1,L2,D);
K_QyQy  =Kern_Ensemble_Sub(Kernel{105},x1.Qy, x2.Qy' ,y1.Qy, y2.Qy', A,L1,L2,D,[],sigma_Qy,BC.Qy);
K_QyMx  =Kern_Ensemble_Sub(Kernel{106},x1.Qy, x2.Mx' ,y1.Qy, y2.Mx', A,L1,L2,D,nu);
K_QyMy  =Kern_Ensemble_Sub(Kernel{107},x1.Qy, x2.My' ,y1.Qy, y2.My', A,L1,L2,D,nu);
K_QyMxy =Kern_Ensemble_Sub(Kernel{108},x1.Qy, x2.Mxy',y1.Qy, y2.Mxy',A,L1,L2,D,nu);

% Row Moment Mx 
K_Mxw   =Kern_Ensemble_Sub(Kernel{109},x1.Mx, x2.w',  y1.Mx, y2.w',  A,L1,L2,D,nu);
K_MxRx  =Kern_Ensemble_Sub(Kernel{110},x1.Mx, x2.Rx', y1.Mx, y2.Rx', A,L1,L2,D,nu);
K_MxRy  =Kern_Ensemble_Sub(Kernel{111},x1.Mx, x2.Ry', y1.Mx, y2.Ry', A,L1,L2,D,nu);
K_MxKx  =Kern_Ensemble_Sub(Kernel{112},x1.Mx, x2.Kx', y1.Mx, y2.Kx', A,L1,L2,D,nu);
K_MxKy  =Kern_Ensemble_Sub(Kernel{113},x1.Mx, x2.Ky', y1.Mx, y2.Ky', A,L1,L2,D,nu);
K_MxKxy =Kern_Ensemble_Sub(Kernel{114},x1.Mx, x2.Kxy',y1.Mx, y2.Kxy',A,L1,L2,D,nu);
K_Mxp   =Kern_Ensemble_Sub(Kernel{115},x1.Mx, x2.p'  ,y1.Mx, y2.p',  A,L1,L2,D,nu);
K_MxQx  =Kern_Ensemble_Sub(Kernel{116},x1.Mx, x2.Qx' ,y1.Mx, y2.Qx', A,L1,L2,D,nu);
K_MxQy  =Kern_Ensemble_Sub(Kernel{117},x1.Mx, x2.Qy' ,y1.Mx, y2.Qy', A,L1,L2,D,nu);
K_MxMx  =Kern_Ensemble_Sub(Kernel{118},x1.Mx, x2.Mx' ,y1.Mx, y2.Mx', A,L1,L2,D,nu,sigma_Mx,BC.Mx);
K_MxMy  =Kern_Ensemble_Sub(Kernel{119},x1.Mx, x2.My' ,y1.Mx, y2.My', A,L1,L2,D,nu);
K_MxMxy =Kern_Ensemble_Sub(Kernel{120},x1.Mx, x2.Mxy',y1.Mx, y2.Mxy',A,L1,L2,D,nu);

% Row Moment My 
K_Myw   =Kern_Ensemble_Sub(Kernel{121},x1.My, x2.w',  y1.My, y2.w',  A,L1,L2,D,nu);
K_MyRx  =Kern_Ensemble_Sub(Kernel{122},x1.My, x2.Rx', y1.My, y2.Rx', A,L1,L2,D,nu);
K_MyRy  =Kern_Ensemble_Sub(Kernel{123},x1.My, x2.Ry', y1.My, y2.Ry', A,L1,L2,D,nu);
K_MyKx  =Kern_Ensemble_Sub(Kernel{124},x1.My, x2.Kx', y1.My, y2.Kx', A,L1,L2,D,nu);
K_MyKy  =Kern_Ensemble_Sub(Kernel{125},x1.My, x2.Ky', y1.My, y2.Ky', A,L1,L2,D,nu);
K_MyKxy =Kern_Ensemble_Sub(Kernel{126},x1.My, x2.Kxy',y1.My, y2.Kxy',A,L1,L2,D,nu);
K_Myp   =Kern_Ensemble_Sub(Kernel{127},x1.My, x2.p'  ,y1.My, y2.p',  A,L1,L2,D,nu);
K_MyQx  =Kern_Ensemble_Sub(Kernel{128},x1.My, x2.Qx' ,y1.My, y2.Qx', A,L1,L2,D,nu);
K_MyQy  =Kern_Ensemble_Sub(Kernel{129},x1.My, x2.Qy' ,y1.My, y2.Qy', A,L1,L2,D,nu);
K_MyMx  =Kern_Ensemble_Sub(Kernel{130},x1.My, x2.Mx' ,y1.My, y2.Mx', A,L1,L2,D,nu);
K_MyMy  =Kern_Ensemble_Sub(Kernel{131},x1.My, x2.My' ,y1.My, y2.My', A,L1,L2,D,nu,sigma_My,BC.My);
K_MyMxy =Kern_Ensemble_Sub(Kernel{132},x1.My, x2.Mxy',y1.My, y2.Mxy',A,L1,L2,D,nu);

% Row Moment Mxy 
K_Mxyw  =Kern_Ensemble_Sub(Kernel{133},x1.Mxy,x2.w',  y1.Mxy,y2.w',  A,L1,L2,D,nu);
K_MxyRx =Kern_Ensemble_Sub(Kernel{134},x1.Mxy,x2.Rx', y1.Mxy,y2.Rx', A,L1,L2,D,nu);
K_MxyRy =Kern_Ensemble_Sub(Kernel{135},x1.Mxy,x2.Ry', y1.Mxy,y2.Ry', A,L1,L2,D,nu);
K_MxyKx =Kern_Ensemble_Sub(Kernel{136},x1.Mxy,x2.Kx', y1.Mxy,y2.Kx', A,L1,L2,D,nu);
K_MxyKy =Kern_Ensemble_Sub(Kernel{137},x1.Mxy,x2.Ky', y1.Mxy,y2.Ky', A,L1,L2,D,nu);
K_MxyKxy=Kern_Ensemble_Sub(Kernel{138},x1.Mxy,x2.Kxy',y1.Mxy,y2.Kxy',A,L1,L2,D,nu);
K_Mxyp  =Kern_Ensemble_Sub(Kernel{139},x1.Mxy,x2.p'  ,y1.Mxy,y2.p',  A,L1,L2,D,nu);
K_MxyQx =Kern_Ensemble_Sub(Kernel{140},x1.Mxy,x2.Qx' ,y1.Mxy,y2.Qx', A,L1,L2,D,nu);
K_MxyQy =Kern_Ensemble_Sub(Kernel{141},x1.Mxy,x2.Qy' ,y1.Mxy,y2.Qy', A,L1,L2,D,nu);
K_MxyMx =Kern_Ensemble_Sub(Kernel{142},x1.Mxy,x2.Mx' ,y1.Mxy,y2.Mx', A,L1,L2,D,nu);
K_MxyMy =Kern_Ensemble_Sub(Kernel{143},x1.Mxy,x2.My' ,y1.Mxy,y2.My', A,L1,L2,D,nu);
K_MxyMxy=Kern_Ensemble_Sub(Kernel{144},x1.Mxy,x2.Mxy',y1.Mxy,y2.Mxy',A,L1,L2,D,nu,sigma_Mxy,BC.Mxy);

Kern  =[K_ww      K_wRx      K_wRy      K_wKx     K_wKy     K_wKxy    K_wp     K_wQx     K_wQy     K_wMx     K_wMy     K_wMxy;...
        K_Rxw     K_RxRx     K_RxRy     K_RxKx    K_RxKy    K_RxKxy   K_Rxp    K_RxQx    K_RxQy    K_RxMx    K_RxMy    K_RxMxy;...
        K_Ryw     K_RyRx     K_RyRy     K_RyKx    K_RyKy    K_RyKxy   K_Ryp    K_RyQx    K_RyQy    K_RyMx    K_RyMy    K_RyMxy;...
        K_Kxw     K_KxRx     K_KxRy     K_KxKx    K_KxKy    K_KxKxy   K_Kxp    K_KxQx    K_KxQy    K_KxMx    K_KxMy    K_KxMxy;...
        K_Kyw     K_KyRx     K_KyRy     K_KyKx    K_KyKy    K_KyKxy   K_Kyp    K_KyQx    K_KyQy    K_KyMx    K_KyMy    K_KyMxy;...
        K_Kxyw    K_KxyRx    K_KxyRy    K_KxyKx   K_KxyKy   K_KxyKxy  K_Kxyp   K_KxyQx   K_KxyQy   K_KxyMx   K_KxyMy   K_KxyMxy;...
        K_pw      K_pRx      K_pRy      K_pKx     K_pKy     K_pKxy    K_pp     K_pQx     K_pQy     K_pMx     K_pMy     K_pMxy;...
        K_Qxw     K_QxRx     K_QxRy     K_QxKx    K_QxKy    K_QxKxy   K_Qxp    K_QxQx    K_QxQy    K_QxMx    K_QxMy    K_QxMxy;...
        K_Qyw     K_QyRx     K_QyRy     K_QyKx    K_QyKy    K_QyKxy   K_Qyp    K_QyQx    K_QyQy    K_QyMx    K_QyMy    K_QyMxy;...
        K_Mxw     K_MxRx     K_MxRy     K_MxKx    K_MxKy    K_MxKxy   K_Mxp    K_MxQx    K_MxQy    K_MxMx    K_MxMy    K_MxMxy;...
        K_Myw     K_MyRx     K_MyRy     K_MyKx    K_MyKy    K_MyKxy   K_Myp    K_MyQx    K_MyQy    K_MyMx    K_MyMy    K_MyMxy;...
        K_Mxyw    K_MxyRx    K_MxyRy    K_MxyKx   K_MxyKy   K_MxyKxy  K_Mxyp   K_MxyQx   K_MxyQy   K_MxyMx   K_MxyMy   K_MxyMxy];

    
end


