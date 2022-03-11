%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% begin with RK2

i1 = 1;
i2 = 2;

Cv_star = zeros(Nz,Nx);
CT_star = zeros(Nz,Nx);
max_adv_nl=zeros(Nt-1,1);
max_adv_l=zeros(Nt-1,1);
vorticity_full=zeros(t(end),Nz/2+1,Nx);
for n = 1:Nt-1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Predictor step
    % --- linear advection

    Cp = squeeze(Cv(i1,:,:)).*inv_laplacian;

    Cwrk = zeros(NzB,NxB);
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv(i1,1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv(i1,Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv(i1,1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv(i1,Nz/2+2:Nz,Nx/2+2:Nx));
    dvortdxB = ifft2(Cwrk,'symmetric')*Ncor;
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv(i1,1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv(i1,Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv(i1,1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv(i1,Nz/2+2:Nz,Nx/2+2:Nx));
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
    
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(CT(i1,1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT(i1,Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT(i1,1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT(i1,Nz/2+2:Nz,Nx/2+2:Nx));
    dthetadxB = ifft2(Cwrk,'symmetric')*Ncor;

    lin_adv = -(U0B + U1B).*dvortdxB - wwindB.*Q0zB - uwindB.*dQ1dxB - wwindB.*dQ1dzB;
    max_adv_l(n)=max(abs(lin_adv),[],'all');
    
    Cwrk = fft2(lin_adv);
    Clin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Clin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Clin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Clin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;
    max_adv_l(n)=max(abs(lin_adv),[],'all');
% --- bouyancy term

    boy = -g*dthetadxB./theta0B;
    Cwrk = fft2(boy);
    Cb(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cb(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cb(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cb(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- eddy viscosity

    Cev = rv*squeeze(Cv(i1,:,:)).*laplacian4;

% --- sponge boundary

    vorticity = ifft2(squeeze(Cv(i1,:,:)),'symmetric');
    Csb = -fft2(sb.*vorticity);

% --- nonlinear advection
    %%add noise to advection
    %usave(jk,n)=max(abs(uwindB),[],'all');
    %uWsave(jk,n)=usave(jk,n)+W(n);
    nonlin_adv = -(uwindB+W(n)).*dvortdxB - wwindB.*dvortdzB;
    max_adv_nl(n)=max(abs(nonlin_adv),[],'all');
    
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- sub-grid scale mixing

    flag_predict = 0;
    CT1 = squeeze( CT(i1,:,:) );
    subgrid_scale_mixing;    
    Dv_n = Dv;
    Dt_n = Dt;
    %Khsave(:,:,i1) = Kh;
    
% --- tendency
    
    Cv_tendency(i1,:,:) = Clin_adv + Cnonlin_adv + Cb + CF + Cev + Csb;

    % integrate one time step

    Cv_star = squeeze(Cv(i1,:,:) + dt*Cv_tendency(i1,:,:)) + dt*Dv_n;
    
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
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(CT(i1,1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT(i1,Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT(i1,1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT(i1,Nz/2+2:Nz,Nx/2+2:Nx));
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

    Ctd = rt*squeeze(CT(i1,:,:)).*laplacian4;

% --- sponge boundary

    theta = ifft2(squeeze(CT(i1,:,:)),'symmetric');
    Csb = -fft2(sb.*theta);

% --- nonlinear advection
    %%add noise to the advection term
    nonlin_adv = -(uwindB+W(n)).*dthetadxB - wwindB.*dthetadzB;
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- tendency

    CT_tendency(i1,:,:) = Clin_adv + Cnonlin_adv + Cb + CH + Ctd + Csb;

    % integrate one time step

    CT_star = squeeze( CT(i1,:,:) + dt*CT_tendency(i1,:,:) ) + dt*Dt_n;
  
    % enforce BC

    theta = ifft2( CT_star,'symmetric');
    theta(1,:) = 0;
    theta(Nz/2+1,:) = 0;
    for k = Nz/2+2:Nz
        theta(k,:) = -theta(Nz-k+2,:);
    end
    CT_star = fft2(theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Corrector step
    % --- linear advection

    Cp = squeeze(Cv_star).*inv_laplacian;

    Cwrk = zeros(NzB,NxB);
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv_star(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv_star(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv_star(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv_star(Nz/2+2:Nz,Nx/2+2:Nx));
    dvortdxB = ifft2(Cwrk,'symmetric')*Ncor;
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(Cv_star(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(Cv_star(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(Cv_star(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(Cv_star(Nz/2+2:Nz,Nx/2+2:Nx));
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
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddx(1:Nz/2+1,1:Nx/2+1).*squeeze(CT_star(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddx(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT_star(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddx(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT_star(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddx(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT_star(Nz/2+2:Nz,Nx/2+2:Nx));
    dthetadxB = ifft2(Cwrk,'symmetric')*Ncor;

    lin_adv = -(U0B + U1B).*dvortdxB - wwindB.*Q0zB - uwindB.*dQ1dxB - wwindB.*dQ1dzB;
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

    Cev = rv*Cv_star.*laplacian4;

% --- sponge boundary

    vorticity = ifft2( Cv_star,'symmetric');
    Csb = -fft2(sb.*vorticity);

% --- nonlinear advection
    %%add noise to the advection
    nonlin_adv = -(uwindB+W(n)).*dvortdxB - wwindB.*dvortdzB;
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- sub-grid scale mixing

    flag_predict = 1;
    CT1 = squeeze( CT_star );
    subgrid_scale_mixing;  
    
% --- tendency

    Cv_tendency(i2,:,:) = Clin_adv + Cnonlin_adv + Cb + CF + Cev + Csb;

    % integrate one time step

    Cv(i2,:,:) = squeeze( Cv(i1,:,:) + dt/2*( Cv_tendency(i2,:,:) + Cv_tendency(i1,:,:) ) )...
                                     + dt * ( (1-alpha_dc)*Dv_n + alpha_dc*Dv );
    
    % enforce BC

    vorticity = ifft2(squeeze(Cv(i2,:,:)),'symmetric');
    vorticity(1,:) = 0;
    vorticity(Nz/2+1,:) = 0;
    for k = Nz/2+2:Nz
        vorticity(k,:) = -vorticity(Nz-k+2,:);
    end
    Cv(i2,:,:) = fft2(vorticity);
    
    % potential temperature equation

    Cwrk = zeros(NzB,NxB);
    Cwrk(1:Nz/2+1,1:Nx/2+1) = ddz(1:Nz/2+1,1:Nx/2+1).*squeeze(CT_star(1:Nz/2+1,1:Nx/2+1));
    Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1) = ddz(Nz/2+2:Nz,1:Nx/2+1).*squeeze(CT_star(Nz/2+2:Nz,1:Nx/2+1));
    Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB) = ddz(1:Nz/2+1,Nx/2+2:Nx).*squeeze(CT_star(1:Nz/2+1,Nx/2+2:Nx));
    Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = ddz(Nz/2+2:Nz,Nx/2+2:Nx).*squeeze(CT_star(Nz/2+2:Nz,Nx/2+2:Nx));
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

    Ctd = rt*CT_star.*laplacian4;

% --- sponge boundary

    theta = ifft2(CT_star,'symmetric');
    Csb = -fft2(sb.*theta);

% --- nonlinear advection
    %%add noise W(n) to the advection
    nonlin_adv = -(uwindB+W(n)).*dthetadxB - wwindB.*dthetadzB;
    Cwrk = fft2(nonlin_adv);
    Cnonlin_adv(1:Nz/2+1,1:Nx/2+1) = Cwrk(1:Nz/2+1,1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,1:Nx/2+1) = Cwrk(NzB-Nz/2+2:NzB,1:Nx/2+1)/Ncor;
    Cnonlin_adv(1:Nz/2+1,Nx/2+2:Nx) = Cwrk(1:Nz/2+1,NxB-Nx/2+2:NxB)/Ncor;
    Cnonlin_adv(Nz/2+2:Nz,Nx/2+2:Nx) = Cwrk(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB)/Ncor;

% --- tendency
    clear uwindB wwindB dthetadxB dthetadzB dvortdxB dvortdzB
    CT_tendency(i2,:,:) = Clin_adv + Cnonlin_adv + Cb + CH + Ctd + Csb;

    % integrate one time step

    CT(i2,:,:) = squeeze( CT(i1,:,:) + dt/2*( CT_tendency(i2,:,:) + CT_tendency(i1,:,:) ) ) ...
                                     + dt * ( (1-alpha_dc)*Dt_n + alpha_dc*Dt ) ;
  
    % enforce BC

    theta = ifft2(squeeze(CT(i2,:,:)),'symmetric');
    theta(1,:) = 0;
    theta(Nz/2+1,:) = 0;
    for k = Nz/2+2:Nz
        theta(k,:) = -theta(Nz-k+2,:);
    end
    CT(i2,:,:) = fft2(theta);
    
    % advance indices
    % advance indices
    %Cv_full(n+1,:,:) = Cv(i2,:,:);
    %CT_full(n+1,:,:) = CT(i2,:,:);
    if(mod(t(n+1),1)==0)
        vorticity_full(t(n+1),:,:) = vorticity(1:Nz/2+1,:);
    end
    theta_full(1,:,:) = theta(1:Nz/2+1,:);
    
    %i2 = i2+1;
    %i1 = i1+1;
    Cv(1,:,:)=Cv(2,:,:);
    CT(1,:,:)=CT(2,:,:);
    
    
    
end
Cv_save=Cv(2,:,:);


