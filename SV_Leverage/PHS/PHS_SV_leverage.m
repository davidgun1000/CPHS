%profile on
%load the dataset
load('dataindustrystockreal3000.mat');
y=y_order_demean_used(6,:);

%load the dataset
%load('leverage_6000.mat');

%initial values for the parameters
phi=0.98;
tau=0.10225;
rho=-0.4592;
mu=-0.4239;
%the number of particles
N=5000;
%the number of MCMC iterations
nloop=15000;
%the length of time series
T=length(y);

%initial values of loglikelihood and log volatilities
[X,W,A,~]=smc_mu_univ_leverage(y,phi,tau,mu,rho,N,T);
ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,y);

s=1;
target_accept=0.25;
D1=2;
V1=0.01*eye(D1);
scale=1;
accept_tau=0;

a_phi=100;
b_phi=1.5;

for i=1:nloop
        i
        rho
        tic
        
        eps=(y(1,1:T-1)).*exp(-ctraj(1,1:T-1)./2);
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
        
        %conditional sequential monte carlo algorithm
        [X,W,A,lik]=csmc_mu_univ_leverage(y,phi,tau,mu,rho,N,T,ctraj);
        prior=logcauchy(tau); %calculate log prior at the current parameters
        post=prior+lik; %calculate log posterior at the current parameters
        jac=log(1/tau); %calculate the log jacobian
        
        theta=[log(tau),atanh(rho)]; 
        R1=mvnrnd(theta,scale.*V1);% adaptive random walk proposals
        tau_star=exp(R1(1,1));
        rho_star=tanh(R1(1,2));
        %SMC algorithm
        [X_star,W_star,A_star,lik_star]=smc_mu_univ_leverage(y,phi,tau_star,mu,rho_star,N,T);
        prior_star=logcauchy(tau_star); %calculate the log prior at the proposed parameters
        post_star=prior_star+lik_star; %calculate the log posterior at the proposed parameters
        jac_star=log(1/tau_star); %calculate the log jacobian at the proposed parameters
        r1=exp(post_star-post+jac-jac_star); %compute the MH acceptance ratio
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
        if i>500
            scale=update_sigma(scale,C1,target_accept,i,2);
            V1=cov(thetasave(1:i,:));
            V1=jitChol(V1);
        end

        ctraj=backward_simulation_SV_leverage_univ(X,W,A,phi,tau,mu,rho,y);

        toc   
        %save the results
         if mod(i,1000)==0
            save('standardMix_N10000_T6000_leverage.mat','Post');         
	     end
        
        Post.phi(i,1)=phi;
        Post.tau(i,1)=tau;
        Post.mu(i,1)=mu;
        Post.rho(i,1)=rho;
        Post.scale(i,1)=scale;
        id=(1:1:300)*10;
        Post.ctraj(i,:)=ctraj(1,id);
        
        

        
        
        
        
end
save('standardMix_N10000_T6000_leverage.mat','Post');
%         [X,W,A,lik]=csmc_mu_univ_leverage(y,phi,tau,mu,rho,N,T,ctraj);