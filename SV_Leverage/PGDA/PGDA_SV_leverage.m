%profile on
%load the dataset
%load('dataindustrystockreal3000.mat');
%y=y_order_demean_used(6,:);

%load the dataset
load('leverage_6000.mat');

%initial values of the parameters
phi=0.98;
tau=0.01;
mu=-0.4820;
rho=-0.4;
phi_trans=logit_inverse(phi);
tau_trans=log(tau);
mu_trans=mu;
rho_trans=atanh(rho);

%the number of particles
N=5000;
nloop=11000; %the number of MCMC iterations
T=length(y); %the length of time series

ctraj=[randn(1,T);phi*ones(1,T);mu*ones(1,T);tau*ones(1,T);rho*ones(1,T)];
num_param=4;

prior_hyper.a0=100;
prior_hyper.b0=1.5;

pseudo_hyper.phi_trans=0.001;
pseudo_hyper.mu_trans=0.001;
pseudo_hyper.tau_trans=0.001;
pseudo_hyper.rho_trans=0.005;

param_trans.phi=phi_trans;
param_trans.mu=mu_trans;
param_trans.tau=tau_trans;
param_trans.rho=rho_trans;

for i=1:nloop
        i
        ctraj(5,1)
        tic
        
        %generating pseudo observations for the parameters
        pseudo_obs.phi_trans=normrnd(phi_trans,sqrt(pseudo_hyper.phi_trans));
        pseudo_obs.mu_trans=normrnd(mu_trans,sqrt(pseudo_hyper.mu_trans));
        pseudo_obs.tau_trans=normrnd(tau_trans,sqrt(pseudo_hyper.tau_trans));
        pseudo_obs.rho_trans=normrnd(rho_trans,sqrt(pseudo_hyper.rho_trans));
        
        %conditional sequential monte carlo used in PGDA algorithm. 
        [X,W,A]=csmc_mu_SV_Fearnhead_leverage(y,pseudo_obs,pseudo_hyper,prior_hyper,param_trans,N,T,ctraj);
        %ancestor tracing algorithm for the PGDA algorithm. 
        [ctraj]=ancestor_trace_Fearnhead_augmentation(X,W,A);
        phi_trans=logit_inverse(ctraj(2,1));
        mu_trans=ctraj(3,1);
        tau_trans=log(ctraj(4,1));
        rho_trans=atanh(ctraj(5,1));
        param_trans.phi=phi_trans;
        param_trans.mu=mu_trans;
        param_trans.tau=tau_trans;
        param_trans.rho=rho_trans;
        Post.phi(i,1)=ctraj(2,1);
        Post.mu(i,1)=ctraj(3,1);
        Post.tau(i,1)=ctraj(4,1);
        Post.rho(i,1)=ctraj(5,1);
        
        Post.phi_trans(i,1)=phi_trans;
        Post.mu_trans(i,1)=mu_trans;
        Post.tau_trans(i,1)=tau_trans;
        Post.rho_trans(i,1)=rho_trans;
        id=(1:1:300)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        if i>=500
           pseudo_hyper.phi_trans=var(Post.phi_trans);
           pseudo_hyper.mu_trans=var(Post.mu_trans);
           pseudo_hyper.tau_trans=var(Post.tau_trans);
           pseudo_hyper.rho_trans=var(Post.rho_trans);  
        end
        %save the results in your directory
        if mod(i,1000)==0
%            save('PG_Fearnhead_leverage_sim_N1000_T6000.mat','Post');
        end

        toc    
        
end
save('PG_Fearnhead_leverage_sim_N1000_T6000.mat','Post');
