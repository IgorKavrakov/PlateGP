function [AP] = Example1_Analytical_SSRecPlate(AP)
% Analytical solution for the plate model in Example 1
x=AP.x(:);y=AP.y(:); q0=AP.q0; a=AP.a; b=AP.b; D=AP.D; nu=AP.nu;
Nx=length(x);Ny=length(y);

Nm=1; Nn=1;
m=1:2:Nm*2;n=1:2:Nn*2;n=n(:);m=m(:);

pmn=16*q0/pi^2./(m*n'); %Load magnitude 
p=sin(x*pi*m'/a)*pmn*sin(n*pi*y'/b);

wmn=pmn./(D*pi^4*(repmat(m,[1 length(n)]).^2/a^2+repmat(n',[length(m) 1]).^2/b^2).^2);%Displacements mag per Fourier
w=sin(x*pi*m'/a)*wmn*sin(n*pi*y'/b); % Displacement

kxmn=(pi^2/a^2)*repmat((m.^2),[1 length(n)]).*wmn;
Kx=sin(x*pi*m'/a)*kxmn*sin(n*pi*y'/b); % Curvature x
kymn=(pi^2/b^2)*repmat((n'.^2),[length(m) 1]).*wmn;
Ky=sin(x*pi*m'/a)*kymn*sin(n*pi*y'/b); % Curvature y
kxymn=(pi^2/a/b)*m*n'.*wmn;
Kxy=-2*cos(x*pi*m'/a)*kxymn*cos(n*pi*y'/b); % Curvature xy

Mx=D.*(Kx+nu.*Ky); % Moment Mx
My=D.*(Ky+nu.*Kx); % Moment My
Mxy=D.*(1-nu).*-Kxy./2; % Moment Mxy

qxmnkX=(pi/a).*repmat((m),[1 length(n)]).*kxmn;
qxmnkY=(pi/a).*repmat((m),[1 length(n)]).*kymn;
Qx=D*(cos(x*pi.*m'/a).*qxmnkX.*sin(n*pi*y'/b)+cos(x*pi.*m'/a).*qxmnkY.*sin(n*pi*y'/b));

qymnkX=(pi/b).*repmat((n'.^2),[length(m) 1]).*kxmn;
qymnkY=(pi/b).*repmat((n'.^2),[length(m) 1]).*kymn;
Qy=D*(sin(x*pi.*m'/a).*qymnkX.*cos(n*pi*y'/b)+sin(x*pi.*m'/a).*qymnkY.*cos(n*pi*y'/b));

% Store output
[AP.x_mesh,AP.y_mesh]=meshgrid(x,y); 
AP.p=p; AP.w=w'; AP.Kx=Kx'; AP.Ky=Ky'; AP.Kxy=Kxy'; AP.Qx=Qx'; AP.Qy=Qy'; AP.Mx=Mx'; AP.My=My'; AP.Mxy=Mxy';

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
end

