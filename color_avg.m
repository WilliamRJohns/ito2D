%Methods to compare with, rk2D_rk2P is the "true" solution
meth1='itoeul'
meth2='rk3rk3'
gamma=15;        %Scaling parameter for the Noise
tic
tF=32;          %length of integration in seconds
Nx=128*2;       %x grid resolution (also needs to be changed in background_state)
Nz = Nx/2;      %zgrid reolution

%arrays for avg convergence over different alpha values
avg_eul=[]; avg_time_eul=[];
avg_ito=[]; avg_time_ito=[];
avg_rk2=[]; avg_time_rk2=[];

%%different colors of noise
%for alpha=[0,1e-5,1e-4,1e-3,1]
for alpha=[0]
    alpha 
    %arrays for different noise realizations
    diff_mat_eul=[]; time_mat_eul=[];
    diff_mat_ito=[]; time_mat_ito=[];
    diff_mat_rk2=[]; time_mat_rk2=[];
    %%number of realizations to average over
    nr=1;
    for ik=1:nr
        %Vectors for final error and running times for different dt resolutions
        rel_diff_eul=[];    time_eul=[];
        rel_diff_ito=[];    time_ito=[];
        rel_diff_rk2=[]; time_rk2=[];
        rel_theta_ito=[];
        rel_theta_rk2=[];
        rel_theta_eul=[];
        
        ik
        dt=.01 % Min dt used for "groud truth" solution
        Nt = ( tF/dt + 1 ); 

        %Generate the white noise and ito correction term
        %gamma=20;
        seed=100*ik;
        rng(seed,'twister')
        t = [0:dt:dt*(Nt-1)]; %time vector
        eta=zeros(Nt-1,2*Nz,2*Nx);
        for zz=1:length(t)-1
            [eta(zz,:,:),Ito_sum]=color_noise_2d_uniform(t(zz),length(t)-1,alpha, seed, Nz, Nx,gamma);
        end
        W=eta/sqrt(dt);
        
        %%Compute high res rk2 solution "ground truth"
        [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,'rk2rk2',W,Ito_sum);
        %time_rk2=[time_rk2 t_end];
        %"Ground Truth" solutions for comparison
        vort_mindt=squeeze(vorticity_full(end,:,:));
        theta_mindt=squeeze(theta_full(end,:,:));
        
%         %run high res Eul+Ito model
%         [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,meth1,W,Ito_sum);
%         time_ito=[time_ito t_end];
%         %Compute final error at tF
%         temp=norm(squeeze(vorticity_full(end,:,:))-vort_mindt)/norm(vort_mindt);
%         rel_diff_ito=[rel_diff_ito temp];
%         temp=norm(squeeze(theta_full(end,:,:))-theta_mindt)/norm(theta_mindt);
%         rel_theta_ito=[rel_theta_ito temp];
%         
%         %run high res Eul model
%         [vorticity_full,theta_full]=run_model_switch(tF,dt,Nx,meth2,W,Ito_sum);
%         time_eul=[time_eul t_end];
%         temp=norm(squeeze(vorticity_full(end,:,:))-vort_mindt)/norm(vort_mindt);
%         rel_diff_eul=[rel_diff_eul temp];
%         temp=norm(squeeze(theta_full(end,:,:))-theta_mindt)/norm(theta_mindt);
%         rel_theta_eul=[rel_theta_eul temp];
       
        % Run the model for increasing values of dt
        steps=[.05,.1,.5,1,2,4,8,16];
        %steps=[.005,.01,.05,.1,.5,1,2,4,8,16];
        for dt=steps
            dt
            %Average the noise over the new dt
            Nt = ( tF/dt + 1 );
            jj=(size(W,1)/(Nt-1));
            clear eta2
            for zz=1:(size(eta,1)-jj)/jj+1
                eta2(zz,:,:)=squeeze(mean(eta(jj*(zz-1)+1:jj*zz,:,:)));
            end
            eta=eta2*sqrt(jj); %Eta is again ~N(0,1)
            W=eta/sqrt(dt);    %W for the new dt
            
            %Run rk2 model 
            [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,'rk2rk2',W,Ito_sum);
            time_rk2=[time_rk2 t_end];
            temp=norm(squeeze(vorticity_full(end,:,:))-vort_mindt)/norm(vort_mindt);
            rel_diff_rk2=[rel_diff_rk2 temp];
            temp=norm(squeeze(theta_full(end,:,:))-theta_mindt)/norm(theta_mindt);
            rel_theta_rk2=[rel_theta_rk2 temp];
            
            %Run Eul+ito model
            [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,meth1,W,Ito_sum);
            time_ito=[time_ito t_end];
            %Compute final error at tF
            temp=norm(squeeze(vorticity_full(end,:,:))-vort_mindt)/norm(vort_mindt);
            rel_diff_ito=[rel_diff_ito temp];
            temp=norm(squeeze(theta_full(end,:,:))-theta_mindt)/norm(theta_mindt);
            rel_theta_ito=[rel_theta_ito temp];
            
            %Run Eul model
            [vorticity_full,theta_full,t_end]=run_model_switch(tF,dt,Nx,meth2,W,Ito_sum);
            time_eul=[time_eul t_end];
            temp=norm(squeeze(vorticity_full(end,:,:))-vort_mindt)/norm(vort_mindt);
            rel_diff_eul=[rel_diff_eul temp];
            temp=norm(squeeze(theta_full(end,:,:))-theta_mindt)/norm(theta_mindt);
            rel_theta_eul=[rel_theta_eul temp];

        end %end dt loop
    %Add errors for this noise realization to diff_matt    
    diff_mat_eul=[diff_mat_eul ;rel_diff_eul];
    save(sprintf('diff_mat_eul_%dw_%d.mat',gamma,alpha),'diff_mat_eul')
    time_mat_eul=[time_mat_eul; time_eul];

    diff_mat_ito=[diff_mat_ito ;rel_diff_ito];
    save(sprintf('diff_mat_ito_%dw_%d.mat',gamma,alpha),'diff_mat_ito')
    time_mat_ito=[time_mat_ito; time_ito];

    diff_mat_rk2=[diff_mat_rk2 ;rel_diff_rk2];
    save(sprintf('diff_mat_rk2_%dw_%d.mat',gamma,alpha),'diff_mat_rk2')
    time_mat_rk2=[time_mat_rk2; time_rk2];
    end %end nr loop (noise realizations
%Compute average errors over all noise realizations    
avg_eul=[avg_eul ; mean(diff_mat_eul,1)];
save(sprintf('avg_eul_%dw.mat',gamma),'avg_eul')
avg_time_eul=[avg_time_eul; mean(time_mat_eul,1)];

avg_ito=[avg_ito ; mean(diff_mat_ito,1)];
save(sprintf('avg_ito_%dw.mat',gamma),'avg_ito')
avg_time_ito=[avg_time_ito; mean(time_mat_ito,1)];

avg_rk2=[avg_rk2 ; mean(diff_mat_rk2,1)];
save(sprintf('avg_rk2_%dw.mat',gamma),'avg_rk2')
avg_time_rk2=[avg_time_rk2; mean(time_mat_rk2,1)];
end
toc
%steps=[.01 steps];

figure()
semilogx(steps,avg_time_rk2(1,:)); hold on;
semilogx(steps,avg_time_eul(1,:));
semilogx(steps,avg_time_ito(1,:));
legend('rk2',meth2,meth1)
title(sprintf('Avg running times alpha=%d'),0);
xlabel('dt')
ylabel('seconds')

figure()
loglog(steps,avg_rk2(1,:),'b-*');hold on
loglog(steps,avg_eul(1,:),'r-x');
loglog(steps,avg_ito(1,:));
legend('rk2',meth2,meth1)
title(sprintf('Covergence alpha=%d,gamma=%d',0,gamma))

figure()
loglog(steps,avg_rk2(2,:));hold on
loglog(steps,avg_eul(2,:));
loglog(steps,avg_ito(2,:));
legend('rk2','eul','ito')
title(sprintf('Covergence alpha=%d',1e-4))

figure()
loglog(steps,avg_rk2(3,:));hold on
loglog(steps,avg_eul(3,:));
loglog(steps,avg_ito(3,:));
legend('rk2','eul','ito')
title(sprintf('Covergence alpha=%d',1e-3))

figure()
loglog(steps,avg_rk2(1,:));hold on
loglog(steps,avg_eul(1,:));
loglog(steps,avg_ito(1,:));
legend('rk2','eul','ito')
title(sprintf('Covergence alpha=%d',1))

figure()
loglog(steps,diff_mat_rk2,'b');hold on
loglog(steps,diff_mat_ito,'r*');
loglog(steps,diff_mat_eul,'y');
legend('rk2','eul','ito')
title(sprintf('Covergence alpha=%d',0))

figure()
loglog(steps,avg_ito)
legend('0','1e-3','1e-4','1')
title('Ito color')

figure()
loglog(steps,avg_eul)
legend('0','1e-3','1e-4','1')
title('Eul color')

figure()
loglog(steps,avg_rk2)
legend('0','1e-3','1e-4','1')
title('rk2 color')

%linear fit curves for convergence order ( asymptotic regime must be
%identified for limits)
rk2fit=polyfit(steps(2:5),avg_rk2(2:5),1);
fprintf('RK2 order %d\n',rk2fit(1))
itofit=polyfit(steps(4:5),avg_ito(4:5),1);
fprintf('Ito order %d\n',itofit(1))

% figure()
% loglog(steps,avg_eul,'-*');hold on
% title(sprintf('Euler Convergence dt=.1 gamma=%d,tf=%d',gamma,tF))
% legend('\alpha=0','\alpha=1e-5','\alpha=1e-4','\alpha=1e-3','\alpha=1')
% xlabel('dt')
% ylabel('rel_err')
% 
% figure()
% loglog(steps,avg_ito,'-o');hold on
% title(sprintf('Ito Convergence dt=.1 gamma=%d,tf=%d',gamma,tF))
% legend('\alpha=0','\alpha=1e-5','\alpha=1e-4','\alpha=1e-3','\alpha=1')
% xlabel('dt')
% ylabel('rel_err')
% 
% figure()
% loglog(steps,avg_rk2,'--');hold on
% title(sprintf('Rk2 Convergence dt=.1 gamma=%d,tf=%d',gamma,tF))
% legend('\alpha=0','\alpha=1e-5','\alpha=1e-4','\alpha=1e-3','\alpha=1')
% xlabel('dt')
% ylabel('rel_err')
% 
% 
% % figure()
% % loglog(steps,diff_mat_rk2(end,:));hold on
% % loglog(steps,diff_mat_eul(end,:));
% % loglog(steps,diff_mat_ito(end,:));
% % legend('rk2','eul','ito')
% % title(sprintf('Covergence alpha=%d',1))
% % 
% figure()
% loglog(steps,diff_mat_rk2(1,:));hold on
% loglog(steps,diff_mat_eul(1,:),'r-*');
% loglog(steps,diff_mat_ito(1,:));
% legend('rk2','eul','ito')
% title(sprintf('Covergence alpha=%d',0))
% % 
% % figure()
% % loglog(steps,diff_mat_rk2(3,:));hold on
% % loglog(steps,diff_mat_eul(3,:));
% % loglog(steps,diff_mat_ito(3,:));
% % legend('rk2','eul','ito')
% % title(sprintf('Covergence alpha=%d',1e-4))
% 
% 
% 
% % figure()
% % loglog(steps,mean(diff_mat_rk2,1),'g-x');hold on
% % loglog(steps,mean(diff_mat_ito,1),'-o')
% % loglog(steps,mean(diff_mat_eul,1),'-*');
% % legend('rk2','ito','eul')
% % title(sprintf('Avg Convergence vort rk2 dt=.1 gamma=%d,tF=%d',gamma,tF))
% % xlabel('dt')
% % ylabel('rel_err')
% % err=zeros(1,Nt);
% % for ii=1:Nt
% %  err(ii)=norm(squeeze(vorticity_eul_full(ii,:,:))-squeeze(vorticity_ito_full(ii,:,:)));
% % end
% % figure()
% % plot(err)
