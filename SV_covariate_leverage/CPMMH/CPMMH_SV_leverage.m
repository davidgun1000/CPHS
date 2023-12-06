%series number 6 for the real data
%load the data
load('leverage_6000_cov.mat');
T=length(y');
%initial values for parameters
phi=phi_true;
tau=tau_true;
mu=mu_true;
rho=rho_true;
%number of covariates
num_beta=50;
beta=beta_true;
%the number of particles 
N=1000;
%the number of MCMC iterations
nloop=15000;

%hyperparameters for phi
a0_phi=100; 
b0_phi=1.5;


u1=randn(N,T);
u1_res=[rand(N,T)];
%initial values of log-likelihood
[X,W,A,lik]=smc_mu_univ_corr_leverage(y,phi,tau,mu,rho,beta,z,N,T,u1,u1_res);
%backward simulation algorithm
ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,beta,z,y);

%calculating the log prior
prior=logcauchy(tau)+log(betapdf((1+phi)/2,a0_phi,b0_phi))+sum(log(mvnpdf(beta',zeros(1,num_beta),eye(num_beta))));
%calculating the log jacobian
jac=log(1/tau)+log((1/phi)+(1/(1-phi)));
%calculating the log posterior
post=prior+lik;

s=1;
%target acceptance probability
target_accept=0.25;

%initial scale and covariance matrix for the random walk proposals
D1=54;
V1=0.001*eye(D1);
scale=0.1;
rho_corr=0.999;
accept=0;
for i=1:nloop
        i
        tic
        A1=rand();
        %proposals for the us
        u1_star=u1*rho_corr+sqrt(1-rho_corr^2)*randn(N,T);
        u1_res_norm=norminv(u1_res);
        u1_res_norm_star=u1_res_norm*rho_corr+sqrt(1-rho_corr^2)*randn(N,T);
        u1_res_star=normcdf(u1_res_norm_star);
        %generating proposals using random walk
        theta=[mu,logit_inverse(phi),log(tau),atanh(rho),beta'];
        R1=mvnrnd(theta,scale.*V1);
        mu_star=R1(1,1);
        phi_star=logitcdf(R1(1,2));
        tau_star=exp(R1(1,3));
        rho_star=tanh(R1(1,4));
        beta_star=R1(1,5:end)';
        %calculating the log prior evaluated at the proposed parameters
        prior_star=logcauchy(tau_star)+log(betapdf((1+phi_star)/2,a0_phi,b0_phi))+sum(log(mvnpdf(beta_star',zeros(1,num_beta),eye(num_beta))));
        [X_star,W_star,A_star,lik_star]=smc_mu_univ_corr_leverage(y,phi_star,tau_star,mu_star,rho_star,beta_star,z,N,T,u1_star,u1_res_star);
        if sum(isnan(lik_star))>0 | sum(isnan(X_star))>0 | sum(isnan(W_star))>0 | sum(isnan(A_star))>0 | rho_star>0.999  
           phi=phi;
           tau=tau;
           mu=mu;
           rho=rho;
           beta=beta;
           post=post;
           jac=jac;
           u1=u1;
           u1_res=u1_res;
           X=X;
           W=W;
           A=A;
           ctraj=ctraj;  
           theta=[mu,logit_inverse(phi),log(tau),atanh(rho),beta'];
           thetasave(i,:)=theta;
        else
        %backward simulation algorithm    
        ctraj_star=backward_simulation_SV_leverage_univ(X_star,W_star,A_star,phi_star,tau_star,mu_star,rho_star,beta_star,z,y);
        %calculating the log posterior
        post_star=prior_star+lik_star;
        %calculating the log jacobian
        jac_star=log(1/tau_star)+log((1/phi_star)+(1/(1-phi_star)));
       
        %calculating the MH ratio
        r1 = exp(post_star-post+jac-jac_star);
        C1 = min(1,r1);   
        if A1<=C1
           phi=phi_star;
           tau=tau_star;
           mu=mu_star;
           rho=rho_star;
           beta=beta_star;
           post=post_star;
           jac=jac_star;
           u1=u1_star;
           u1_res=u1_res_star;
           accept=accept+1;
           ctraj=ctraj_star;
        end
        theta=[mu,logit_inverse(phi),log(tau),atanh(rho),beta'];
        thetasave(i,:)=theta;
        
        end
        
        if i>50
             scale=update_sigma(scale,C1,target_accept,i,D1);
             V1=cov(thetasave(1:i,:));
             V1=jitChol(V1);
        end
        
        toc
        
        
         %save the results
         if mod(i,1000)==0
            %save('/srv/scratch/z3512791/univ_SV/corrPMMH_N20_T3000_leverage.mat','Post');
            %save('/short/jz21/dg2271/corr_mix/corrPMMH_N50_T3000_realbeta.mat','Post');
            %save('corrPMMH_N10_T3000_realbeta.mat','Post');
            save('/scratch/jz21/dg2271/corrPMMH_sim_N1000_T6000_cov.mat','Post');
         end
        
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        Post.beta(i,:)=beta;
        id=(1:1:T/10)*10;
        Post.ctraj(i,:)=ctraj(1,id);
          
end
save('/scratch/jz21/dg2271/corrPMMH_sim_N1000_T6000_cov.mat','Post');
% save('series3partperiod_leverage1pc.mat','Post');% 2 demaen


%eta_mean_store=mean(Post.eta(5000:end,:))';
%eta_median_store=prctile(Post.eta(5000:end,:),50)';

%end

