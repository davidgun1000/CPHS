%code for CPMMH


%load the real data
%load('dataindustrystockreal3000.mat');
%y=y_order_demean_used(6,:);

%load the simulated data
load('leverage_6000.mat');
T=length(y); %the length of time series
phi=0.98; %initial values
tau=0.10;
mu=-0.42;
rho=-0.45;
N=100; %the number of particles
nloop=15000; %the number of MCMC iterations

%hyperparameters for the prior of phi
a0_phi=100; 
b0_phi=1.5;

%initial values for the log-likelihood
u1=randn(N,T);
u1_res=[rand(N,T)];
[X,W,A,lik]=smc_mu_univ_corr_leverage(y,phi,tau,mu,rho,N,T,u1,u1_res);
%generating the initial log-trajectory 
ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,y);

%calculating the prior density
prior=logcauchy(tau)+log(betapdf((1+phi)/2,a0_phi,b0_phi));
%calculating the jacobian
jac=log(1/tau)+log((1/phi)+(1/(1-phi)));
%calculating the log-posterior
post=prior+lik;

s=1;
%target MH acceptance probability
target_accept=0.25;
%initial values for the covariance matrix
D1=4;
V1=0.01*eye(D1);
scale=1;
rho_corr=0.999;
accept=0;

for i=1:nloop
        i
        tic
        A1=rand();
        %update the random numbers
        u1_star=u1*rho_corr+sqrt(1-rho_corr^2)*randn(N,T);
        u1_res_norm=norminv(u1_res);
        u1_res_norm_star=u1_res_norm*rho_corr+sqrt(1-rho_corr^2)*randn(N,T);
        u1_res_star=normcdf(u1_res_norm_star);
        %generate the proposal for the parameters
        theta=[mu,logit_inverse(phi),log(tau),atanh(rho)];
        R1=mvnrnd(theta,scale.*V1);
        mu_star=R1(1,1);
        phi_star=logitcdf(R1(1,2));
        tau_star=exp(R1(1,3));
        rho_star=tanh(R1(1,4));
        %calculating the log prior
        prior_star=logcauchy(tau_star)+log(betapdf((1+phi_star)/2,a0_phi,b0_phi));
        %calculating the log-likelihood using particle filter
        [X_star,W_star,A_star,lik_star]=smc_mu_univ_corr_leverage(y,phi_star,tau_star,mu_star,rho_star,N,T,u1_star,u1_res_star);
        
        if sum(isnan(lik_star))>0 | sum(isnan(X_star))>0 | sum(isnan(W_star))>0 | sum(isnan(A_star))>0 | rho_star>0.999  
           phi=phi;
           tau=tau;
           mu=mu;
           rho=rho;
           post=post;
           jac=jac;
           u1=u1;
           u1_res=u1_res;
           X=X;
           W=W;
           A=A;
           ctraj=ctraj;  
           theta=[mu,logit_inverse(phi),log(tau),atanh(rho)];
           thetasave(i,:)=theta;
        else
        %backward simulation algorithm     
        ctraj_star=backward_simulation_SV_leverage_univ(X_star,W_star,A_star,phi_star,tau_star,mu_star,rho_star,y);
        %calculating the log posterior evaluated at the proposed parameters
        post_star=prior_star+lik_star;
        %calculating the log-jacobian evaluated at the proposed parameters
        jac_star=log(1/tau_star)+log((1/phi_star)+(1/(1-phi_star)));
       
        %calculating the MH acceptance ratio
        r1 = exp(post_star-post+jac-jac_star);
        C1 = min(1,r1);   
        if A1<=C1
           phi=phi_star;
           tau=tau_star;
           mu=mu_star;
           rho=rho_star;
           post=post_star;
           jac=jac_star;
           u1=u1_star;
           u1_res=u1_res_star;
           accept=accept+1;
           ctraj=ctraj_star;
        end
        theta=[mu,logit_inverse(phi),log(tau),atanh(rho)];
        thetasave(i,:)=theta;
        
        end
        
        %updating the covariance matrix of the proposals 
        if i>50
             scale=update_sigma(scale,C1,target_accept,i,5);
             V1=cov(thetasave(1:i,:));
             V1=jitChol(V1);
        end
        
        
        
        
         %save the results
         if mod(i,1000)==0
            save('corrPMMH_sim_N100_T6000.mat','Post');
         end
        
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        id=(1:1:T/10)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        toc  
end
%save('/scratch/jz21/dg2271/corrPMMH_sim_N100_T6000.mat','Post');

% save('series3partperiod_leverage1pc.mat','Post');% 2 demaen


%eta_mean_store=mean(Post.eta(5000:end,:))';
%eta_median_store=prctile(Post.eta(5000:end,:),50)';

%end

