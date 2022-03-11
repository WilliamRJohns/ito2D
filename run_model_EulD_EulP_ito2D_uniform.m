%%Euler Dynamics

Cv_star = zeros(Nz,Nx);
CT_star = zeros(Nz,Nx);
Cv = squeeze(Cv_0(1,:,:));
CT = squeeze(CT_0(1,:,:));
Itov=zeros(Nz,Nx);
Itot=zeros(Nz,Nx);
CIto=zeros(Nz,Nx);
vorticity_full=zeros(t(end),Nz/2+1,Nx);
max_adv_nl=zeros(Nt-1,1);
max_adv_l=zeros(Nt-1,1);
for n = 1:Nt-1
    % --- linear advection

    Cp = Cv.*inv_laplacian;

    Cwrk = zeros(NzB,NxB);
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv(Nz/2+2:Nz,Nx/2+2:Nx));
    dvortdxB = ifft2(Cwrk,'symmetric')*Ncor;
    
    %Compute dvort^2/dx^2 for ito correction
    CwrkI = zeros(NzB,NxB);
    CwrkI(1:Nz/2+1,1:Nx/2+1)= ddx(1:Nz/2+1,1:Nx/2+1).*Cwrk(1:Nz/2+1,1:Nx/2+1);
    CwrkI(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1);
    CwrkI(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB);
    CwrkI(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB);
    dvort2dx2B = ifft2(CwrkI,'symmetric')*Ncor;
    %test(n,:,:)=dvort2dx2B;
    
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv(Nz/2+2:Nz,Nx/2+2:Nx));
    dvortdzB = ifft2(Cwrk,'symmetric')*Ncor;
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*Cp(1:Nz/2+1,1:Nx/2+1);
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*Cp(Nz/2+2:Nz,1:Nx/2+1);
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*Cp(1:Nz/2+1,Nx/2+2:Nx);
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*Cp(Nz/2+2:Nz,Nx/2+2:Nx);
    uwindB = ifft2(Cwrk,'symmetric')*Ncor;
    Cwrk(1:Nz/2+1,1:Nx/2+1) = -ddx(1:Nz/2+1,1:Nx/2+1).*Cp(1:Nz/2+1,1:Nx/2+1);
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = -ddx(Nz/2+2:Nz,1:Nx/2+1).*Cp(Nz/2+2:Nz,1:Nx/2+1);
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = -ddx(1:Nz/2+1,Nx/2+2:Nx).*Cp(1:Nz/2+1,Nx/2+2:Nx);
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = -ddx(Nz/2+2:Nz,Nx/2+2:Nx).*Cp(Nz/2+2:Nz,Nx/2+2:Nx);
    wwindB = ifft2(Cwrk,'symmetric')*Ncor;
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(CT(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT(Nz/2+2:Nz,Nx/2+2:Nx));
    dthetadxB = ifft2(Cwrk,'symmetric')*Ncor;
    
    %Compute dtheta^2/dx^2 for Ito correction
    CwrkI = zeros(NzB,NxB);
    CwrkI(1:Nz/2+1,1:Nx/2+1)= ddx(1:Nz/2+1,1:Nx/2+1).*Cwrk(1:Nz/2+1,1:Nx/2+1);
    CwrkI(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1);
    CwrkI(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB);
    CwrkI(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB);
    dtheta2dx2B = ifft2(CwrkI,'symmetric')*Ncor;
    %test(1,:,:)=dtheta2dx2B;
 
    lin_adv = -(U0B + U1B).*dvortdxB - wwindB.*Q0zB - uwindB.*dQ1dxB - wwindB.*dQ1dzB;
    max_adv_l(n)=max(abs(lin_adv),[],'all');
    
    Cwrk = fft2(lin_adv);
    Clin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Clin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Clin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Clin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- bouyancy term

    boy = -g*dthetadxB./theta0B;
    Cwrk = fft2(boy);
    Cb(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cb(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cb(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cb(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- eddy viscosity

    Cev = rv*squeeze(Cv(:,:)).*laplacian4;

% --- sponge boundary

    vorticity = ifft2(squeeze(Cv(:,:)),'symmetric');
    Csb = -fft2(sb.*vorticity);

% --- nonlinear advection

    nonlin_adv = -(uwindB+W(n)).*dvortdxB - wwindB.*dvortdzB;
    max_adv_nl(n)=max(abs(nonlin_adv),[],'all');
    
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- sub-grid scale mixing

    flag_predict = 0;
    CT1 = CT;
    subgrid_scale_mixing;    
    Dv_n = Dv;
    Dt_n = Dt;
    %Km_full(n,:,:) = Km;
    %CKm_full(n,:,:) = fft2(Km);
    %Dv_full(n,:,:) = Dv;
    %Dt_full(n,:,:) = Dt;
% --- tendency

    Cv_tendency = Clin_adv + Cnonlin_adv + Cb + CF + Cev + Csb;

    % integrate one time step

    Cv_star = Cv + dt*Cv_tendency + dt*Dv_n;
    
    %Ito correction
    Cwrk = fft2(dvort2dx2B);
    CIto(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    CIto(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    CIto(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    CIto(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;
    Itov=Ito_sum*CIto;
    Cv_star=Cv_star+dt*Itov;
   
    % enforce BC
    vorticity = ifft2( Cv_star,'symmetric');
    vorticity(1,:) = 0;
    vorticity(Nz/2+1,:) = 0;
    for k = Nz/2+2:Nz
        vorticity(k,:) = -vorticity(Nz-k+2,:);
    end
    Cv_star = fft2(vorticity);
    
    % potential temperature equation

    Cwrk = zeros(NzB,NxB);
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(CT(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT(Nz/2+2:Nz,Nx/2+2:Nx));
    dthetadzB = ifft2(Cwrk,'symmetric')*Ncor;

% --- linear advection

    lin_adv = -(U0B + U1B).*dthetadxB - W1B.*dthetadzB - uwindB.*dtheta1dxB - wwindB.*dtheta1dzB;
    Cwrk = fft2(lin_adv);
    Clin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Clin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Clin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Clin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- bouyancy term

    boy = -wwindB.*theta0zB;
    Cwrk = fft2(boy);
    Cb(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cb(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cb(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cb(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- thermal diffusion

    Ctd = rt*CT.*laplacian4;

% --- sponge boundary

    theta = ifft2(CT,'symmetric');
    Csb = -fft2(sb.*theta);

% --- nonlinear advection

    nonlin_adv = -(uwindB+W(n)).*dthetadxB - wwindB.*dthetadzB;
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- tendency

    CT_tendency = Clin_adv + Cnonlin_adv + Cb + CH + Ctd + Csb;

    % integrate one time step

    CT_star = CT + dt*CT_tendency  + dt*Dt_n;
    
    %Ito correction
    
    Cwrk = fft2(dtheta2dx2B);
    CIto(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    CIto(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    CIto(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    CIto(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;
    Itot=Ito_sum*CIto;
    CT_star=CT_star+dt*Itot;
    % enforce BC

    theta = ifft2( CT_star,'symmetric');
    theta(1,:) = 0;
    theta(Nz/2+1,:) = 0;
    for k = Nz/2+2:Nz
        theta(k,:) = -theta(Nz-k+2,:);
    end
    CT_star = fft2(theta);
    
    
    % advance indices
    if(mod(t(n+1),1)==0)
        vorticity_full(t(n+1),:,:) = vorticity(1:Nz/2+1,:);
    end
    theta_full(n+1,1:Nz/2,:) = theta(1:Nz/2,:);
    Cv=Cv_star;
    CT=CT_star;
    %Cv_full(n+1,:,:) = Cv_star;
    %CT_full(n+1,:,:) = CT_star;
    
end


