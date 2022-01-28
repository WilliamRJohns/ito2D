function [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,model,W,Ito_sum)
%

%% Parameters
var_factor = 1e-6;
ito_factor = 0;

%Nx = 128*2; % x grid resolution
T0 = 305; % surface temperature
N02 = 1e-4; % bouyancy
sb_factor = 0; % turn on sponge boundaries
Amp0 = 0.01; % initial amplitude of perturbations
kmix = 0.21;% sub-grid mixing parameter
alpha_dc = 0.5; % de-centering
mu = 0.0025; % diffusion parameter

%% dimensions
V = 10; % meters/second
g = 9.8; % meters/second^2
L = 1000; % meters
Fr2 = V^2/(g*L); % square Froude number

%% Initial State and preallocate arrays
background_state
initialization0
initialization    
clear theta_full_0 vorticity_full_0   CT_full_0

%% Run the selected Model
t_start=cputime;
switch model
    case 'rk2rk2'       
        run_model_RK2D_RK2P_sto_adv_u_2D;
    case 'eulrk2'
        run_model_EulD_RK2P_sto_adv_2D
    case 'itork2'
        run_model_EulD_RK2P_ito2D
    case 'euleul'
        run_model_EulD_EulP_sto_adv_2D
    case 'itoeul'
        run_model_EulD_EulP_ito2D
    case 'rk3rk3'
        run_model_RK3D_RK3P_sto_adv_u_2D
    otherwise
        warning('Specifeid Model was not found')
        break
end %end model switch
t_end=cputime-t_start;
end