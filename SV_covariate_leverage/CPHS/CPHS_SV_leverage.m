%profile on
%load the dataset
load('leverage_6000_cov.mat');
%initial values of the parameters
phi=0.97;
tau=0.05;
mu=-0.2;
rho=-0.2;
%the number of covariates
num_beta=50;
%the number of particles
N=100;
%the number of MCMC iterations
nloop=15000;
%the length of time series
T=length(y');

%initial values for the log-volatilities
ctraj=randn(1,T);
s=1;

%target MH acceptance
target_accept=0.25;
%initial covariance matrix for the random walk proposals
D1=2;
V1=0.001*eye(D1);
%initial scale for the random walk proposals
scale=1;
accept_tau=0;

%initial values for the beta
beta=0*ones(num_beta,1);
a_phi=100;
b_phi=1.5;

for i=1:nloop
        i   
        tic
        
        %sampling beta in PG step
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
        %sampling mu in PG step
        mu_var=(tau(1,1)*(1-rho(1,1)^2))/((1-rho(1,1)^2)*(1-phi(1,1)^2)+(T-1)*(1-phi(1,1))^2);
        mu_mean=(mu_var/(tau(1,1)*(1-rho(1,1)^2)))*((1-rho(1,1)^2)*(1-phi(1,1)^2)*(ctraj(1,1))+(1-phi(1,1)).*sum(ctraj(1,2:T)-phi(1,1).*(ctraj(1,1:T-1))-rho(1,1).*sqrt(tau(1,1)).*eps));
        mu(1,1)=normrnd(mu_mean,sqrt(mu_var));
         
        %sampling phi in PG step
        phi_var=(tau(1,1)*(1-rho(1,1)^2))/(sum((ctraj(1,1:T-1)-mu(1,1)).^2)-((ctraj(1,1)-mu(1,1))^2)*(1-rho(1,1)^2));
        phi_mean=sum((ctraj(1,2:T)-mu(1,1)).*(ctraj(1,1:T-1)-mu(1,1))-rho(1,1).*sqrt(tau(1,1)).*(ctraj(1,1:T-1)-mu(1,1)).*eps)./(sum((ctraj(1,1:T-1)-mu(1,1)).^2)-((ctraj(1,1)-mu(1,1))^2)*(1-rho(1,1)^2));
        if phi_mean>0.999
           phi_mean=0.999;
        end
        phi_star=normt_rnd(phi_mean,phi_var,0,0.99999);
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
        
        %conditional sequential Monte Carlo algorithm
        [X,W,A,lik,u1,u1_res_rand]=csmc_mu_univ_corr_leverage(y,phi,tau,mu,rho,beta,z,N,T,ctraj);
         prior=logcauchy(tau); %calculating the log prior
         
         post=prior+lik; %calculating the log posterior
         jac=log(1/tau); %calculating the log jacobian
         
         theta=[log(tau),atanh(rho)];
         R1=mvnrnd(theta,scale.*V1); %random walk proposals
         tau_star=exp(R1(1,1));
         rho_star=tanh(R1(1,2));
         if rho_star>0.999 
            rho_star=0.999;              
         end
         if rho_star<-0.999
            rho_star=-0.999; 
         end
         
         %sequential Monte Carlo algorithm
         [X_star,W_star,A_star,lik_star]=smc_mu_univ_corr_leverage(y,phi,tau_star,mu,rho_star,beta,z,N,T,u1,u1_res_rand);
         prior_star=logcauchy(tau_star); %calculating the log prior evaluated at the proposed parameters
         
         post_star=prior_star+lik_star;  %calculating the log posterior evaluated at the proposed parameters
         jac_star=log(1/tau_star);  %calculating the log jacobian evaluated at the proposed parameters
         
         r1=exp(post_star-post+jac-jac_star); %calculating MH acceptance ratio
         C1=min(1,r1);
         A1=rand();
         if A1<=C1
            tau=tau_star;           
            rho=rho_star;
            accept_tau=accept_tau+1;
            X=X_star;
            W=W_star;
            A=A_star;
         end
         
         thetasave(i,:)=[log(tau),atanh(rho)];
         if i>250
            scale=update_sigma(scale,C1,target_accept,i,2);
            V1=cov(thetasave(1:i,:));
            V1=jitChol(V1);
         end
         %backward simulation algorithm
         ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,beta,z,y);
%        save the results
         toc  
          if mod(i,1000)==0
             save('/scratch/jz21/dg2271/corrmixPMMH_sim_N100_T6000_cov.mat','Post');
          end
%         
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        Post.beta(i,:) = beta';

        Post.scale(i,1)=scale;
        id=(1:1:T/10)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        
        
        
        
        
        
        
end
save('/scratch/jz21/dg2271/corrPMMH_sim_N100_T6000.mat','Post');
 
%         tic
%         z_reshape = reshape(z,num_beta,1,T);
%         z_reshape_trans = multitransp(z_reshape);
%         ctraj_reshape_trans=(reshape(ctraj,1,1,T));
%         z_z_transpose=multiprod(z_reshape,z_reshape_trans);
%         var_temp1 = sum(z_z_transpose./exp(ctraj_reshape_trans),3);
%         var_temp2 = sum(((rho^2).*tau.*multiprod(exp(-ctraj_reshape_trans(:,:,1:T-1)),z_z_transpose(:,:,1:T-1)))./(tau*(1-rho^2)),3);
%         var_beta = inv(var_temp1+var_temp2+eye(num_beta));
%         
%         y_reshape_trans=reshape(y,1,1,T);
%         mean_temp1 = sum(multiprod(z_reshape,y_reshape_trans)./exp(ctraj_reshape_trans),3);
%         
%         mean_temp2 = (sum((rho*sqrt(tau).*exp(-ctraj_reshape_trans(:,:,1:T-1)./2).*(multiprod(-(ctraj_reshape_trans(1,1,2:T)-mu),z_reshape_trans(:,:,1:T-1))+...
%             phi.*multiprod((ctraj_reshape_trans(1,1,2:T)-mu),z_reshape_trans(:,:,1:T-1)) + ...
%             rho*sqrt(tau).*exp(-ctraj_reshape_trans(:,:,1:T-1)./2).*y_reshape_trans(:,:,1:T-1).*z_reshape_trans(:,:,1:T-1)))./(tau*(1-rho^2)),3))';
%         mean_beta = var_beta*(mean_temp1+mean_temp2);
%         beta = (mvnrnd(mean_beta,var_beta))';
%         toc
