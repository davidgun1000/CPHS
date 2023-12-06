%profile on
%load the data
load('leverage_6000_cov.mat');
%initial values for the parameters
phi=phi_true;
tau=tau_true;
mu=mu_true;
rho=rho_true;
beta=beta_true;
%the number of covariates
num_beta=50;
%the number of particles
N=1000;
%the number of MCMC iterations
nloop=15000;
%the length of time series
T=length(y');

%initial values for the log-volatilities
ctraj=ctraj_true;
s=1;
%the target acceptance rate
target_accept=0.25;
%initial values for the scale and covariance matrix in the random walk
%proposals
D1=2;
V1=0.001*eye(D1);
scale=1;
accept_tau=0;

%prior hyperparameters
a_phi=100;
b_phi=1.5;

for i=1:nloop
        i
        tic
        
        %sampling beta using PG step
        var_temp1=0;
        var_temp2=0;
        mean_temp1=0;
        mean_temp2=0;
        for t=1:T
            temp1=(z(:,t)*z(:,t)')./exp(ctraj(1,t));
            var_temp1=var_temp1+temp1;
        
            temp1_mean = (z(:,t)*y(1,t))./exp(ctraj(1,t));
            mean_temp1 = mean_temp1 + temp1_mean;
        
        end
        for t=2:T
            temp2=(rho.^2.*tau.*exp(-ctraj(1,t-1)).*(z(:,t-1)*z(:,t-1)'))./(tau*(1-rho^2));
            var_temp2=var_temp2+temp2;
            
            temp2_mean = (rho*sqrt(tau)*exp(-ctraj(1,t-1)./2).*(-(ctraj(1,t)-mu).*z(:,t-1) + phi.*(ctraj(1,t)-mu).*z(:,t-1) + rho*sqrt(tau)*exp(-ctraj(1,t-1)/2)*y(1,t-1)*z(:,t-1)))./(tau*(1-rho^2));
            mean_temp2 = mean_temp2 + temp2_mean;
            
        end
        var_beta = inv(var_temp1+var_temp2+eye(num_beta));
        mean_beta = var_beta*(mean_temp1 + mean_temp2);
        beta = (mvnrnd(mean_beta,var_beta))';
        
        
        eps=(y(1,1:T-1)-beta'*z(:,1:T-1)).*exp(-ctraj(1,1:T-1)./2);
        %sampling mu using PG step
        mu_var=(tau(1,1)*(1-rho(1,1)^2))/((1-rho(1,1)^2)*(1-phi(1,1)^2)+(T-1)*(1-phi(1,1))^2);
        mu_mean=(mu_var/(tau(1,1)*(1-rho(1,1)^2)))*((1-rho(1,1)^2)*(1-phi(1,1)^2)*(ctraj(1,1))+(1-phi(1,1)).*sum(ctraj(1,2:T)-phi(1,1).*(ctraj(1,1:T-1))-rho(1,1).*sqrt(tau(1,1)).*eps));
        mu(1,1)=normrnd(mu_mean,sqrt(mu_var));
        %sampling phi using PG step
        phi_var=(tau(1,1)*(1-rho(1,1)^2))/(sum((ctraj(1,1:T-1)-mu(1,1)).^2)-((ctraj(1,1)-mu(1,1))^2)*(1-rho(1,1)^2));
        phi_mean=sum((ctraj(1,2:T)-mu(1,1)).*(ctraj(1,1:T-1)-mu(1,1))-rho(1,1).*sqrt(tau(1,1)).*(ctraj(1,1:T-1)-mu(1,1)).*eps)./(sum((ctraj(1,1:T-1)-mu(1,1)).^2)-((ctraj(1,1)-mu(1,1))^2)*(1-rho(1,1)^2));
        if phi_mean>0.999
           phi_mean=0.999;
        end
        phi_star=normt_rnd(phi_mean,phi_var,0,0.999);
        prior=log(betapdf((1+phi(1,1))/2,a_phi,b_phi));
        prior_star=log(betapdf((1+phi_star)/2,a_phi,b_phi));
        num_phi=log(sqrt(1-phi_star^2));
        den_phi=log(sqrt(1-phi(1,1)^2));
        r2=exp(num_phi-den_phi+prior_star-prior);
        C2=min(1,r2);
        A2=rand();
        if A2<=C2
           phi(1,1)=phi_star;
        end
        %sampling tau and rho using PG step
        prior=logcauchy(tau(1,1));
        lik=-((T-1)/2)*log(1-rho(1,1)^2)-0.5.*(1./(tau(1,1)*(1-rho(1,1)^2))).*sum(((ctraj(1,2:T)-mu(1,1))-phi(1,1).*(ctraj(1,1:T-1)-mu(1,1))-rho(1,1).*sqrt(tau(1,1)).*eps).^2)-...
            (T/2)*log(tau(1,1))-0.5.*((1-phi(1,1)^2)/tau(1,1)).*((ctraj(1,1)-mu(1,1)).^2);
        post=prior+lik;
        jac=log(1/tau(1,1));
        theta=[log(tau(1,1)),atanh(rho(1,1))];
        theta_star=mvnrnd(theta,scale.*V1(:,:,1));
        tau_star=exp(theta_star(1,1));
        rho_star=tanh(theta_star(1,2));
        prior_star=logcauchy(tau_star);
        lik_star=-((T-1)/2)*log(1-rho_star^2)-0.5.*(1./(tau_star*(1-rho_star^2))).*sum(((ctraj(1,2:T)-mu(1,1))-phi(1,1).*(ctraj(1,1:T-1)-mu(1,1))-rho_star.*sqrt(tau_star).*eps).^2)-...
            (T/2)*log(tau_star)-0.5.*((1-phi(1,1)^2)/tau_star).*((ctraj(1,1)-mu(s,1)).^2);
        post_star=prior_star+lik_star;
        jac_star=log(1/tau_star);
        r2=exp(post_star-post+jac-jac_star);
        C2=min(1,r2);
        A2=rand();
        if A2<=C2
           rho(1,1)=rho_star;
           tau(1,1)=tau_star;
        end
        thetasave{1,1}(i,:)=[log(tau(1,1)),atanh(rho(1,1))];
        if i>50
           V1(:,:,1)=cov(thetasave{1,1}(1:i,:));
           V1(:,:,1)=jitChol(V1(:,:,1));
           scale(1,1)=update_sigma(scale(1,1),C2,target_accept,i,2); 
        end
        [X,W,A,lik]=csmc_mu_univ_leverage(y,phi,tau,mu,rho,beta,z,N,T,ctraj);
        ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,beta,z,y);
        %save the results
         toc
         if mod(i,1000)==0
            save('/scratch/jz21/dg2271/PG_N1000_T6000_leverage_cov.mat','Post');   
	     end
        
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        Post.scale(i,1)=scale;
        Post.beta(i,:)=beta;
        id=(1:1:600)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        
        
end
save('/scratch/jz21/dg2271/PG_N1000_T6000_leverage_cov.mat','Post');   


