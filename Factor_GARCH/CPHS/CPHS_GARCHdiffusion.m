%PMMH pf
%parpool(28)
%load the data
load('dataindustrystockreal3000.mat');

%the length of time series
T=3001;
%the dimension of stock returns
dim_y=26;
y=y_order_demean_used(1:dim_y,:);

%the number of factors
num_fact=4;
%the number of MCMC iterations
burn=5000;
nit=10000;
iter=burn+nit;
%the number of particles
N=100;
delta=1;
%the number of Euler approximation
M=10;
D1=3;

B0_deep=1;

%prior hyperparameters
hp_sig2=1;
v0=10;
s0=1;
v0_sig2=10;
s0_sig2=10;
a_phi=100;
b_phi=1.5;

%initial values for the parameters
a=0.05*ones(dim_y,1);
tau=0.05*ones(dim_y,1);
mu=zeros(dim_y,1);
B(:,1)=rand(dim_y,1);
B(:,2)=[0;rand(dim_y-1,1)];
B(:,3)=[0;0;rand(dim_y-2,1)];
B(:,4)=[0;0;0;rand(dim_y-3,1)];
tau_factor=1./random('gam',v0/2,2/s0,num_fact,1);
phi_factor=0.98*ones(num_fact,1);

%initial scale and covariance matrix for the random walk proposals

for jj=1:dim_y
    scale_idio(jj,1)=1;
    V1(:,:,jj)=0.0001*eye(D1);
end

D2=1;
for jj=1:num_fact
    scale_factor(jj,1)=1;
    V2(:,:,jj)=0.0001*eye(D2);
end

fact_score=randn(num_fact,T);
%initial values of the log volatilities
parfor s=1:dim_y
   [X,W,A,~]=smc_diffusion_GARCH_factor_SV_mean_approx(y(s,:),T,M,N,delta,a(s,1),tau(s,1),mu(s,1),B(s,:),fact_score,s,dim_y,0,0);
   ctraj_idio(s,:)=ancestor_trace(X,W,A);
end
parfor s=dim_y+1:dim_y+num_fact
   [X,W,A,~]=smc_diffusion_GARCH_factor_SV_mean_approx(0,T,M,N,0,0,0,0,0,fact_score,s,dim_y,phi_factor(s-dim_y,1),tau_factor(s-dim_y,1));
   ctraj_factor(s-dim_y,:)=ancestor_trace(X,W,A);
end


accept1=0;
accept2=zeros(dim_y,1);

target_accept=0.20;
for i=1:iter
    i
    %sampling the factor loading matrix
    tic
    length_ctraj_idio=length(ctraj_idio(1,:));
    id_idio=1:10:length_ctraj_idio;
    ctraj_idio_used=ctraj_idio(:,id_idio);
    for s=1:dim_y
        s_bar=min(s,num_fact);
        F=fact_score(1:s_bar,:)';
        arg_Vi_inv=1./(exp(ctraj_idio_used(s,:)));
        Vi_inv=diag(arg_Vi_inv);
        B_var=inv((F'*Vi_inv*F)+eye(s_bar));
        [B_var]=jitChol(B_var);
        chol_B=chol(B_var);
        B_mean=B_var*F'*((Vi_inv*y(s,:)'));
        B(s,1:s_bar)=(mvnrnd(B_mean',chol_B'*chol_B))';        
    end
           
     %deep interweaving of Kaster et al (2017)
     k1=dim_y;
     for s=1:num_fact
       ind=find(abs(B(:,s))==max(abs(B(:,s))));
       max_B=(B(s,s));
       B_star=B(s:dim_y,s)./max_B;
       fact_star=max_B.*fact_score(s,:);
       ctraj_trans=ctraj_factor(s,:)+2*log(abs(max_B));
       mu_fact_old=log(max_B^2);
       mean_prop=(sum(ctraj_trans(1,2:T-1))+((ctraj_trans(1,T)-phi_factor(s,1)*ctraj_trans(1,1))./(1-phi_factor(s,1))))./(T+1/B0_deep);
       var_prop=(tau_factor(s,1)/((1-phi_factor(s,1))^2))/(T+1/B0_deep);
       mu_fact_new=normrnd(mean_prop,sqrt(var_prop));
       A2=rand();
       comp1_old=logmvnpdf(B_star',zeros(1,k1),hp_sig2*exp(-mu_fact_old)*eye(k1));
       comp2_old=logmvnpdf(ctraj_trans(1,1),mu_fact_old,tau_factor(s,1)/(1-phi_factor(s,1)^2));       
       comp3_old=log(exp(mu_fact_old/2-exp(mu_fact_old)/(2*hp_sig2)));
       comp4_old=logmvnpdf(mu_fact_old,0,B0_deep*tau_factor(s,1)/((1-phi_factor(s,1))^2));
       comp1_new=logmvnpdf(B_star',zeros(1,k1),hp_sig2*exp(-mu_fact_new)*eye(k1));
       comp2_new=logmvnpdf(ctraj_trans(1,1),mu_fact_new,tau_factor(s,1)/(1-phi_factor(s,1)^2));
       comp3_new=log(exp(mu_fact_new/2-exp(mu_fact_new)/(2*hp_sig2)));
       comp4_new=logmvnpdf(mu_fact_new,0,B0_deep*tau_factor(s,1)/((1-phi_factor(s,1))^2));
       R2=exp(comp1_new+comp2_new+comp3_new-comp1_old-comp2_old-comp3_old+comp4_old-comp4_new);
       C2=min(1,R2);
       if A2<=C2
          lam1=exp(mu_fact_new/2);
          accept1=accept1+1;
          
       else
           lam1=max_B;
       end     
       B(:,s)=(lam1/max_B).*B(:,s);
       fact_score(s,:)=(max_B/lam1).*fact_score(s,:);
       ctraj_factor(s,:)=ctraj_factor(s,:)+2*log(abs(max_B/lam1));
       k1=k1-1;
     end

    %sampling the latent factor ft for t=1,...,T.
    for t=1:T
        ctraj_idio_temp(:,t)=ctraj_idio_used(1:dim_y,t);
        ctraj_fact_temp(:,t)=ctraj_factor(1:num_fact,t);
        Vt_inv=diag(1./(exp(ctraj_idio_temp(:,t))));
        Dt_inv=diag(exp(-ctraj_fact_temp(:,t)));
        var_ft=inv((B'*Vt_inv*B)+Dt_inv);
        [var_ft]=jitChol(var_ft);
        chol_var_ft=chol(var_ft);
        mean_ft=var_ft*B'*(Vt_inv*y(:,t));
        fact_score(:,t)=mvnrnd(mean_ft,chol_var_ft'*chol_var_ft);      
    end
    

     %sampling the GARCH parameters and generating latent volatilities
     parfor s=1:dim_y
         [X,W,A,lik,u1,u1_res]=csmc_diffusion_GARCH_factor_SV_mean_approx_corr(y(s,:),T,M,N,delta,a(s,1),tau(s,1),mu(s,1),B(s,:),fact_score,ctraj_idio(s,:),s,dim_y,0,0);
         prior=log_IG_PDF_used(a(s,1),v0/2,s0/2)+log_IG_PDF_used(tau(s,1),v0/2,s0/2);
         
         post=prior+lik;
         jac=log(1/a(s,1))+log(1/tau(s,1));
         theta=[log(a(s,1)),log(tau(s,1)),mu(s,1)];
         R1=mvnrnd(theta,scale_idio(s,1).*V1(:,:,s));
         %R1=mvnrnd(theta,scale_idio(s,1).*covmat_idio(:,:,s));
         
         a_star=exp(R1(1,1));
         tau_star=exp(R1(1,2));
         mu_star=R1(1,3);
         [Xstar,Wstar,Astar,lik_star]=smc_diffusion_GARCH_factor_SV_mean_approx_corr(y(s,:),T,M,N,delta,a_star,tau_star,mu_star,B(s,:),fact_score,s,dim_y,0,0,u1,u1_res);
         if sum(isnan(Wstar))==0 & sum(isnan(W))==0
         %ctraj(s,:)=ancestor_trace(X,W,A);
         %ctraj_star=ancestor_trace(Xstar,Wstar,Astar);
         prior_star=log_IG_PDF_used(a_star,v0/2,s0/2)+log_IG_PDF_used(tau_star,v0/2,s0/2);
         
         post_star=prior_star+lik_star;
         jac_star=log(1/a_star)+log(1/tau_star);
         r1 = exp(post_star-post+jac-jac_star);
         C1 = min(1,r1);
         A1=rand();
         if A1<=C1
           a(s,1)=a_star;
           tau(s,1)=tau_star;
           mu(s,1)=mu_star;
           %ctraj(s,:)=ctraj_star;
           X=Xstar;
           W=Wstar;
           A=Astar;
           accept2(s,1)=accept2(s,1)+1;
         end
         thetasave(i,:,s)=[log(a(s,1)),log(tau(s,1)),mu(s,1)];        
         if i>50
            scale_idio(s,1)=update_sigma(scale_idio(s,1),C1,target_accept,i,3); 
         end
         
         ctraj_idio(s,:)=ancestor_trace(X,W,A);
         
         else
            a(s,1)=a(s,1);
            tau(s,1)=tau(s,1);
            mu(s,1)=mu(s,1);
            ctraj_idio(s,:)=ctraj_idio(s,:);
            thetasave(i,:,s)=[log(a(s,1)),log(tau(s,1)),mu(s,1)];
            scale_idio(s,1)=scale_idio(s,1);
            %accept_nan(s,1)=accept_nan(s,1)+1;
         end 
     end
     for s=1:dim_y
     if i>100
        V1(:,:,s)=cov(thetasave(1:i,:,s));
     end
     end
    
    %sampling phi in PG step
    
    for s=1:num_fact
       phi_var=1./(sum((ctraj_factor(s,2:T-1).^2)./tau_factor(s,1)));
       phi_mean=phi_var.*(sum((ctraj_factor(s,2:T).*ctraj_factor(s,1:T-1))./tau_factor(s,1)));
       if phi_mean>1
          phi_mean=1;
       end
       phi_star=normt_rnd(phi_mean,phi_var,0,0.99999);
       num_phi=log(sqrt(1-phi_star^2));
       den_phi=log(sqrt(1-phi_factor(s,1)^2));
       prior=log(betapdf((1+phi_factor(s,1))/2,a_phi,b_phi));
       prior_star=log(betapdf((1+phi_star)/2,a_phi,b_phi));
       r2=exp(num_phi-den_phi+prior_star-prior);
       C2=min(1,r2);
       A2=rand();
       if A2<=C2
          phi_factor(s,1)=phi_star;
       end
    end
    

    
    parfor s=dim_y+1:dim_y+num_fact
           [X,W,A,lik,u1,u1_res]=csmc_diffusion_GARCH_factor_SV_mean_approx_corr(0,T,M,N,delta,0,0,0,0,fact_score,ctraj_factor(s-dim_y,:),s,dim_y,phi_factor(s-dim_y,1),tau_factor(s-dim_y,1)); 
           prior_tau_factor = log_IG_PDF_used(tau_factor(s-dim_y,1),v0/2,s0/2);
           post=prior_tau_factor+lik;
           jac=log(1/tau_factor(s-dim_y,1));
           theta=log(tau_factor(s-dim_y,1));
           R1=mvnrnd(theta,scale_factor(s-dim_y,1).*V2(:,:,s-dim_y));
           tau_factor_star=exp(R1(1,1));
           [X_star,W_star,A_star,lik_star]=smc_diffusion_GARCH_factor_SV_mean_approx_corr(0,T,M,N,delta,0,0,0,0,fact_score,s,dim_y,phi_factor(s-dim_y,1),tau_factor_star,u1,u1_res);
           if sum(isnan(W_star))==0 & sum(isnan(W))==0
           prior_tau_factor_star=log_IG_PDF_used(tau_factor_star,v0/2,s0/2);
           post_star=prior_tau_factor_star+lik_star;
           jac_star = log(1/tau_factor_star);
           r1 = exp(post_star-post+jac-jac_star);
            C1 = min(1,r1);
            A1=rand();
            if A1<=C1
               tau_factor(s-dim_y,1)=tau_factor_star;
               X=X_star;
               W=W_star;
               A=A_star;
            end
            ctraj_factor(s-dim_y,:)=ancestor_trace(X,W,A); 
            thetasave2(i,:,s-dim_y)=log(tau_factor(s-dim_y,1)); 
            
            if i>50
               scale_factor(s-dim_y,1)=update_sigma(scale_factor(s-dim_y,1),C1,target_accept,i,1); 
            end
            
           end
  
    end
    
     for s=dim_y+1:dim_y+num_fact
        if i>250
           V2(:,:,s-dim_y)=cov(thetasave2(1:i,:,s-dim_y));
        end
     end


          
    Post.B1(i,:)=B(:,1)';
    Post.B2(i,:)=B(:,2)';
    Post.B3(i,:)=B(:,3)';
    Post.B4(i,:)=B(:,4)';
    Post.tau(i,:)=tau';
    Post.a(i,:)=a';
    Post.phi(i,:)=phi_factor;
    Post.tau_factor(i,:)=tau_factor;
    Post.mu(i,:)=mu';
    Post.scale(i,:)=scale_idio';
    toc
    %save the results
     if i==250 | i==500 | i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==10000 | i==12000 | i==14000
        save('/scratch/rl75/dg2271/corrmix_GARCH_PHS.mat','Post');
     end
    
end
save('/scratch/rl75/dg2271/corrmix_GARCH_PHS.mat','Post');


