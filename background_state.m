%%
clear U0 U0z Q0Z N2 Int_N2 theta0 theta0z psi1_0 theta1_0
% time domain
Nt = ( tF/dt + 1 ); 
t = [0:dt:dt*(Nt-1)];

% spatial domain

% grid
Nz = Nx/2;
NxB = 2*Nx;
NzB = 2*Nz;
Ncor = 4;

Ndim = Nx*( Nz/2+1 );

xS = 0;
xF = 4*pi*L;
dx = (xF - xS)/(Nx - 0);
x = [xS:dx:xF-dx];

zS = 0;
zF = 2*L;
% NOTE: the z direction includes a "mirror region"
% solution is in z = 0,1 (L?- William)
dz = (zF - zS)/(Nz - 0);
z = zS:dz:zF-dz;
zp = zS:dz:zF/2;

% eddy viscosity and thermal diffusion
% fourth-order hyper-diffusion
kmax = 2*pi/(xF-xS)*Nx/2;
rv = -mu/(kmax^4)/dt; % 4th order diffusion requires a minus sign ...
rt = rv;

CFL = ( V*kmax )*dt;

%disp([' Froude = ',num2str(sqrt(Fr2))])
%disp([' CFL = ',num2str(CFL)])
%disp([' Diffusion parameter = ',num2str(rv)])
%disp([' x domain length (meters) = ', num2str(xF)])
%disp([' x grid resolution (meters) = ', num2str(dx)])
%disp([' z domain length (meters) = ', num2str(zF/2)])
%disp([' z grid resolution (meters) = ', num2str(dz)])
%disp([' length of integration (seconds) = ', num2str(tF)])

%%
% create Fourier factors

deriv_factors;


%%
% background reference state

b = 8;

for k = 1:Nz/2+1
    
    for j = 1:Nx
    
        U0(k,j)  = V*(1 + tanh(b*(z(k)/L - 0.5)))/2;
        U0z(k,j) = V*(b/L)*sech(b*(z(k)/L - 0.5))^2/2; 
        Q0z(k,j) = V/(L^2)*(-2*b^2*sech(b*(z(k)/L - 0.5))^2*tanh(b*(z(k)/L - 0.5)))/2;

% Int_N2 is the vertical integral of N2         
        
        N2(k,j) = N02;
        Int_N2(k,j) = N02*z(k)/L;
        theta0(k,j) = T0*exp(Int_N2(k,j));
        theta0z(k,j) = N2(k,j)*theta0(k,j)/g;
       
    end
end

% make background fields symmetric

for k = Nz/2+2:Nz
    U0(k,:)      = U0(Nz-k+2,:);
    U0z(k,:)     = U0z(Nz-k+2,:);
    Q0z(k,:)     = Q0z(Nz-k+2,:);
    theta0(k,:)  = theta0(Nz-k+2,:);
    theta0z(k,:) = theta0z(Nz-k+2,:);
end

% produce fields on anti-alias grid

C = fft2(U0);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
U0B = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(Q0z);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
Q0zB = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(theta0);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
theta0B = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(theta0z);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
theta0zB = ifft2(Cwrk,'symmetric')*Ncor;

% create heat source

CH = zeros(Nz,Nx);
%H = zeros(Nz,Nx);
CF = zeros(Nz,Nx);
%F = zeros(Nz,Nx);

% comment this in if there is no forced state

for k = 1:Nz
    for j = 1:Nx
       
        U1(k,j) = 0;
        W1(k,j) = 0;
        dQ1dx(k,j) = 0;
        dQ1dz(k,j) = 0;
        dtheta1dx(k,j) = 0;
        dtheta1dz(k,j) = 0;
       
    end
end

% produce fields on anti-alias grid

C = fft2(U1);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
U1B = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(W1);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
W1B = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(dQ1dx);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
dQ1dxB = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(dQ1dz);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
dQ1dzB = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(dtheta1dx);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
dtheta1dxB = ifft2(Cwrk,'symmetric')*Ncor;
C = fft2(dtheta1dz);
Cwrk = zeros(NzB,NxB);
Cwrk(1:Nz/2+1,1:Nx/2+1) = (C(1:Nz/2+1,1:Nx/2+1));
Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = (C(Nz/2+2:Nz,1:Nx/2+1));
Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = (C(1:Nz/2+1,Nx/2+2:Nx));
Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = (C(Nz/2+2:Nz,Nx/2+2:Nx));
dtheta1dzB = ifft2(Cwrk,'symmetric')*Ncor;

%%
% define sponge boundary

sb = zeros(Nz,Nx);

for k = 1:Nz
    for j = 1:Nx
        sb(k,j) = sb_factor*( exp(-2/L^2*(x(j)-xF).^2) + exp(-2/L^2*(x(j)-xS).^2) );
    end
end

Nmodes = Nx/2;