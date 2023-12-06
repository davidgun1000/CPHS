%profile on
%load the dataset
%load('dataindustrystockreal3000.mat');
%y=y_order_demean_used(6,:);

%load the dataset
load('leverage_6000.mat');
%initial values for the parameters
phi=0.98;
tau=0.01;
mu=-0.4820;
rho=-0.1;
%the number of particles
N=1000;
%the number of MCMC iterations
nloop=15000;
%the length of time series
T=length(y);
%initial values for the log trajectory obtained using SMC. 
[X,W,A,~]=smc_mu_univ_leverage(y,phi,tau,mu,rho,N,T);
%backward simulation algorithm
ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,y);
s=1;
%target MH ratio
target_accept=0.25;

%initial covariance matrix for random walk proposals.
D1=2;
V1=0.001*eye(D1);
%intial scale of the covariance matrix of the random walk proposals. 
scale=1;
accept_tau=0;

%hyperparameters for phi
a_phi=100;
b_phi=1.5;

for i=1:nloop
        i
        rho
        tic
        
        eps=(y(1,1:T-1)).*exp(-ctraj(1,1:T-1)./2);
        
        %sampling mu conditioning on the states
        mu_var=(tau(1,1)*(1-rho(1,1)^2))/((1-rho(1,1)^2)*(1-phi(1,1)^2)+(T-1)*(1-phi(1,1))^2);
        mu_mean=(mu_var/(tau(1,1)*(1-rho(1,1)^2)))*((1-rho(1,1)^2)*(1-phi(1,1)^2)*(ctraj(1,1))+(1-phi(1,1)).*sum(ctraj(1,2:T)-phi(1,1).*(ctraj(1,1:T-1))-rho(1,1).*sqrt(tau(1,1)).*eps));
        mu(1,1)=normrnd(mu_mean,sqrt(mu_var));
        
        %sampling phi conditioning on the states
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
        %sampling tau2 and rho conditioning on the states
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
        %conditional sequential Monte Carlo algorithm
        [X,W,A,lik]=csmc_mu_univ_leverage(y,phi,tau,mu,rho,N,T,ctraj);
        %Backward simulation algorithm 
        ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,y);

        toc   
         if mod(i,1000)==0
         %if i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==8000 | i==9000 | i==10000
            %save('/short/jz21/dg2271/corr_mix/univ_SV_leverage/corrPMMHPG_N20_T3000_realleverage.mat','Post');
            save('PG_N1000_T3000_leverage.mat','Post');
            
	     end
        
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        Post.scale(i,1)=scale;
        id=(1:1:300)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        
        
end
save('PG_N1000_T3000_leverage.mat','Post');
