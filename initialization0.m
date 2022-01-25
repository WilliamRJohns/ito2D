%%
% initialize arrays

clear Cv CT Cv_tendency CT_tendency
Cv_0 = zeros(4,Nz,Nx);
CT_0 = zeros(4,Nz,Nx);
Cv_tendency_0 = zeros(4,Nz,Nx);
CT_tendency_0 = zeros(4,Nz,Nx);

theta1_0 = zeros(Nz,Nx);
psi1_0 = zeros(Nz,Nx);

%Cv_full_0 = zeros(Nt,Nz,Nx);
%CT_full_0 = zeros(Nt,Nz,Nx);

theta_full_0 = zeros(Nt,Nz/2+1,Nx);
vorticity_full_0 = zeros(Nt,Nz/2+1,Nx);
%A_full_0 = zeros(Nt-1,Nz,Nx);
%B_full_0 = zeros(Nt-1,Nz,Nx);
%Ri_full_0 = zeros(Nt-1,Nz/2+1,Nx);
%Km_full_0 = zeros(Nt-1,Nz,Nx);
%CKm_full_0 = Km_full_0;
%Dv_full_0 = zeros(Nt-1,Nz,Nx);
%Dt_full_0 = zeros(Nt-1,Nz,Nx);

%Avec = randn(Nmodes,1);
%Bvec = rand(Nmodes,1);
% Trying to remove strange rel_err_time problem
% A,B saved with myseed=1234 as in original Hodyss model 
load('A.mat');
load('B.mat');

for k = 1:Nz/2+1
    for j = 1:Nx

        for jj = 1:Nmodes
         psi1_0(k,j) = psi1_0(k,j) + sqrt( Amp0 )*(1+0.2*Avec(jj))/(jj^3) ...
                                *sin(2*pi/(zF-zS)*z(k))*sin(jj*2*pi/(xF-xS)*(x(j)-xF*Bvec(jj)));
         theta1_0(k,j) = theta1_0(k,j) + sqrt( Amp0 )*(1+0.2*Avec(jj))/(jj^3)...
                                *sin(2*pi/(zF-zS)*z(k))*cos(jj*2*pi/(xF-xS)*(x(j)-xF*Bvec(jj)));
        end

    end
end


% --- enforce BC

psi1_0(1,:)        = 0;
psi1_0(Nz/2+1,:)   = 0;
theta1_0(1,:)      = 0;
theta1_0(Nz/2+1,:) = 0;
for k = Nz/2+2:Nz
    psi1_0(k,:)   = -psi1_0(Nz-k+2,:);
    theta1_0(k,:) = -theta1_0(Nz-k+2,:);
end

Cp_0 = fft2(psi1_0);
Cv_0(1,:,:) = laplacian.*Cp_0;
CT_0(1,:,:) = fft2(theta1_0);

Cv_full_0(1,:,:) = Cv_0(1,:,:);
CT_full_0(1,:,:) = CT_0(1,:,:);

theta_full_0(1,:,:) = theta1_0(1:Nz/2+1,:);
vorticity = ifft2(squeeze(Cv_0(1,:,:)),'symmetric');
vorticity_full_0(1,:,:) = vorticity(1:Nz/2+1,:);