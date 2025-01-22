function [AP] = Example2_Analytical_FixedRecPlate(AP)
% Analytical solution for plate model in example 2
x=AP.x(:);y=AP.y(:); q0=AP.q0; a=AP.a; b=AP.b; D=AP.D; nu=AP.nu;
Nx=length(x);Ny=length(y);

% Fixed samples
Nm=200; Nn=200;
m=1:1:Nm;n=1:1:Nn;n=n(:);m=m(:);
Nmn=Nm*Nn;

if 1<0
    K2=repmat(diag(2*m'.^4),[Nn Nn]);
    K1=zeros(Nmn,Nmn);K3=zeros(Nmn,Nmn);
    for i=1:Nn
    K1(1+(i-1)*Nm:i*Nm,1+(i-1)*Nm:i*Nm)=diag((m.^2+(a/b).^2*i.^2).^2);
    K3(1+(i-1)*Nm:i*Nm,1+(i-1)*Nm:i*Nm)=2*(a/b*i).^4;
    end
    K=K1+K2+K3;
    wmn=K\(q0*a^4/(4*pi^4*D).*ones(Nmn,1));
    clear K K1 K2 K3;
    save('Example2_FixedRecPlate/wmn','wmn');
else
    % Load pre-saved values for speed
    load('Example2_FixedRecPlate/wmn');
end

wmn=reshape(wmn,[Nm,Nn]);
w=(1-cos(2*x/a*pi*m'))*wmn*(1-cos(2*n*pi*y'/b));% Displacement x

rxmn=2*(pi/a)*repmat((m),[1 length(m)]).*wmn;
Rx=(sin(2*x/a*pi*m'))*rxmn*(1-cos(2*n*pi*y'/b));% Rotation x

rymn=2*(pi/b)*repmat((n),[1 length(n)]).*wmn;
Ry=(1-cos(2*x/a*pi*m'))*rymn*(sin(2*n*pi*y'/b));% Rotation y

kxmn=4*(pi^2/a^2)*repmat((m.^2),[1 length(n)]).*wmn;
Kx=-cos(2*x/a*pi*m'/a)*kxmn*(1-cos(2*n*pi*y'/b));% Curvature x

kymn=4*(pi^2/b^2)*repmat((n'.^2),[length(m) 1]).*wmn;
Ky=-(1-cos(2*x/a*pi*m'))*kymn*cos(2*n*pi*y'/b);% Curvature y

kxymn=4*(pi^2/a/b)*m*n'.*wmn;
Kxy=-2*sin(2*x/a*pi*m')*kxymn*sin(2*n*pi*y'/b); % Curvature xy

p=ones(length(x),length(y))*q0;                 % Uniform load

Mx=D.*(Kx+nu.*Ky); %Moment Mx
My=D.*(Ky+nu.*Kx); %Moment My
Mxy=D.*(1-nu).*Kxy./2; %Moment Mxy

% Store output
[AP.x_mesh,AP.y_mesh]=meshgrid(x,y);
AP.p=p; AP.w=w'; AP.Rx=Rx'; AP.Ry=Ry'; AP.Kx=Kx'; AP.Ky=Ky'; AP.Kxy=Kxy'; AP.Qx=[]; AP.Qy=[]; AP.Mx=Mx'; AP.My=My'; AP.Mxy=Mxy';

% Contaminate data additionally with noise -> Store noise separately (depends on Boundary conditions if it is used)
if isfield(AP,'SNR')
    Names=fieldnames(AP.SNR);
    for i=1:numel(Names)
        if AP.SNR.(Names{i})~=0
            Sig=AP.(Names{i})(:); 
            N_Sig=length(Sig);
            Noise=randn([N_Sig,1])*std(Sig)/AP.SNR.(Names{i});
            AP.Noise.(Names{i})=reshape(Noise,[Nx,Ny]);
        end
    end
end

