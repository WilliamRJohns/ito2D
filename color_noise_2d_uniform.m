function [nt,Ito_sum] = color_noise_2d_uniform(t,Nmode,alpha,seed, Nz, Nx,gamma)
%Nmode=1/dt;
Nmode = Nx/2;
T = 1;
j = 1:Nmode;
wf_list = 2*pi/T*j'; %w_m's
NzB=2*Nz;
NxB=2*Nx;

%rng(seed);
b0 = randn/sqrt(2); % b_0/sqrt(2)
%nt_Im = randn(Nz,Nx)/sqrt(2);
Ito_sum=.5;
gt=b0;
for ii = 1:Nmode
    %rng(seed+ii);
    am = randn;
    bm = randn;
    w_alpha = exp(-alpha*wf_list(ii,1)^2); %C(w_m)
    gt = gt + w_alpha * cos(wf_list(ii,1)*t)*bm +...
    w_alpha * sin(wf_list(ii,1)*t)*am;
    Ito_sum=Ito_sum+w_alpha^2;
   
end


% nt=zeros(NzB,NxB);
% nt(1:Nz/2+1,1:Nx/2+1) = gt(1:Nz/2+1,1:Nx/2+1);
% nt(NzB-Nz/2+2:NzB,1:Nx/2+1) = gt(Nz/2+2:Nz,1:Nx/2+1);
% nt(1:Nz/2+1,NxB-Nx/2+2:NxB) = gt(1:Nz/2+1,Nx/2+2:Nx);
% nt(NzB-Nz/2+2:NzB,NxB-Nx/2+2:NxB) = gt(Nz/2+2:Nz,Nx/2+2:Nx);
nt=gt/sqrt(Nmode)*gamma;
Ito_sum=gamma^2*Ito_sum/2/Nmode;
