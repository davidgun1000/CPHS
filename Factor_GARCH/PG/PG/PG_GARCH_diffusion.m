%PMMH pf
%parpool(28)
%load the data
load('dataindustrystockreal3000.mat');
T=3001;%the length of time series
dim_y=26; %the dimension of the stock returns
num_fact=4; %the number of factors
y=y_order_demean_used(1:dim_y,:); 

%initial values for the parameters
a=0.05*ones(dim_y,1);
tau=0.05*ones(dim_y,1);
mu=zeros(dim_y,1);
B(:,1)=rand(dim_y,1);
B(:,2)=[0;rand(dim_y-1,1)];
B(:,3)=[0;0;rand(dim_y-2,1)];
B(:,4)=[0;0;0;rand(dim_y-3,1)];

%prior hyperparameters
v0=10;
s0=1;
a_phi=100;
b_phi=1.5;
tau_factor=1./random('gam',v0/2,2/s0,num_fact,1);
phi_factor=0.98*ones(num_fact,1);

%the number of MCMC iterations
burn=5000;
nit=10000;
iter=burn+nit;
N=1000; %the number of particles

%initial scale and covariance matrix for random walk proposals.
D1=3;
for jj=1:dim_y
    scale_idio(jj,1)=1;
    V1(:,:,jj)=0.001*eye(D1);
end

D2=1;
for jj=1:num_fact
    V2(:,:,jj)=0.001*eye(D2);
end
%initial values for the factors
fact_score=randn(num_fact,T);

B0_deep=1;
hp_sig2=1;


v0_sig2=10;
s0_sig2=10;

accept1=0;
M=10; %the number of Euler steps
delta=1;
h=delta/M;

%initial values for the log-volatilities
parfor s=1:dim_y
   [X,W,A,~]=smc_diffusion_GARCH_factor_SV_mean_approx(y(s,:),T,M,N,delta,a(s,1),tau(s,1),mu(s,1),B(s,:),fact_score,s,dim_y,phi_factor,tau_factor);
    ctraj_idio(s,:)=ancestor_trace(X,W,A);
end

parfor s=1:num_fact
    [X,W,A,lik]=smc_diffusion_GARCH_factor_SV_mean_approx(0,T,M,N,delta,0,0,0,0,fact_score,s+dim_y,dim_y,phi_factor(s,1),tau_factor(s,1)); 
    ctraj_factor(s,:)=ancestor_trace(X,W,A);
end
accept_idio=zeros(dim_y,1);
accept_nan=zeros(dim_y,1);
target_accept=0.20;
for i=1:iter
    i
    %sampling the factor loading matrix in PG step
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
    %deep interweaving algorithm of Kastner et al,
    
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

    %sampling the latent factor ft for all t=1,...,T.
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
        
    %sampling GARCH parameters using pseudo marginal steps.
    for s=1:dim_y
        
         A1=rand();
         prior=log_IG_PDF_used(a(s,1),v0/2,s0/2)+log_IG_PDF_used(tau(s,1),v0/2,s0/2)+log(unifpdf(mu(s,1),-10,10));
         mean_h=ctraj_idio(s,1:end-1)+h.*(a(s,1).*(mu(s,1)-exp(ctraj_idio(s,1:end-1))).*exp(-ctraj_idio(s,1:end-1))-(tau(s,1)/2));
         var_h=h*tau(s,1);
         lik=sum(log(normpdf(ctraj_idio(s,2:end)',mean_h',sqrt(var_h))));
         jac=log(1/a(s,1))+log(1/tau(s,1));
         post=prior+lik;
         theta=[log(a(s,1)),log(tau(s,1)),mu(s,1)];
         R1=mvnrnd(theta,scale_idio(s,1).*V1(:,:,s));
         a_star=exp(R1(1,1));
         tau_star=exp(R1(1,2));
         mu_star=R1(1,3);
         prior_star=log_IG_PDF_used(a_star,v0/2,s0/2)+log_IG_PDF_used(tau_star,v0/2,s0/2)+log(unifpdf(mu_star,-10,10));
         %mean_h_star=ctraj_idio(s,1:end-1)+a_star*(mu_star-ctraj_idio(s,1:end-1))*h;
         mean_h_star=ctraj_idio(s,1:end-1)+h.*(a_star.*(mu_star-exp(ctraj_idio(s,1:end-1))).*exp(-ctraj_idio(s,1:end-1))-(tau_star/2));
         var_h_star=h*tau_star;
         lik_star=sum(log(normpdf(ctraj_idio(s,2:end)',mean_h_star',sqrt(var_h_star))));
         post_star=prior_star+lik_star;
         jac_star=log(1/a_star)+log(1/tau_star);
         r1 = exp(post_star-post+jac-jac_star);
         C1 = min(1,r1);
         if A1<=C1
             a(s,1)=a_star;
             mu(s,1)=mu_star;
             tau(s,1)=tau_star;
         end
         thetasave(i,:,s)=[log(a(s,1)),log(tau(s,1)),mu(s,1)];
         if i>50
             V1(:,:,s)=cov(thetasave(1:i,:,s));
             scale_idio(s,1)=update_sigma(scale_idio(s,1),C1,target_accept,i,3); 
         end
    end
    
   
    

    
     parfor s=1:dim_y
         [X,W,A,~]=csmc_diffusion_GARCH_factor_SV_mean_approx(y(s,:),T,M,N,delta,a(s,1),tau(s,1),mu(s,1),B(s,:),fact_score,ctraj_idio(s,:),s,dim_y,0,0);
         
         if sum(isnan(W))==0
            ctraj_idio(s,:)=ancestor_trace(X,W,A);         
         else
            ctraj_idio(s,:)=ctraj_idio(s,:);
            accept_nan(s,1)=accept_nan(s,1)+1;
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
       
        v1_factor=v0+T;
        s1_factor=s0+(1-phi_factor(s,1)^2)*ctraj_factor(s,1)+sum((ctraj_factor(s,2:T)-phi_factor(s,1).*ctraj_factor(s,1:T-1)).^2);
        tau_factor(s,1)=1./random('gam',v1_factor/2,2./s1_factor);
    
    end
    
    
    
    parfor s=1:num_fact
           [X,W,A,lik]=csmc_diffusion_GARCH_factor_SV_mean_approx(0,T,M,N,delta,0,0,0,0,fact_score,ctraj_factor(s,:),s+dim_y,dim_y,phi_factor(s,1),tau_factor(s,1)); 
           ctraj_factor(s,:)=ancestor_trace(X,W,A);
           
    end
    

      %keyboard  
    Post.B1(i,:)=B(:,1)';
    Post.B2(i,:)=B(:,2)';
    Post.B3(i,:)=B(:,3)';
    Post.B4(i,:)=B(:,4)';

    Post.tau(i,:)=tau';
    Post.a(i,:)=a';
    Post.phi(i,:)=phi_factor;
    Post.tau_factor(i,:)=tau_factor;
    Post.mu(i,:)=mu';
    toc
    
    %save results.
     if i==250 | i==500 | i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==10000 | i==12000 | i==14000
        save(['/scratch/rl75/dg2271/PG_GARCH',numtostr(N),'.mat'],'Post');
     end
    
end
save(['/scratch/rl75/dg2271/PG_GARCH',numtostr(N),'.mat'],'Post');

%     for s=1:dim_y
%         v1=v0+(T-1);
%         s1=s0+((2*a(s,1))/(1-exp(-2*a(s,1))))*sum((ctraj(s,2:end)-exp(-a(s,1))*ctraj(s,1:T-1)).^2);
%         tau(s,1)=1./random('gam',v1/2,2/s1);
% 
%         A1=rand();
%         prior=log(unifpdf(a(s,1),0,100));
%         mean_h=exp(-a(s,1))*ctraj(s,1:T-1);
%         var_h=tau(s,1)*((1-exp(-2*a(s,1)))/(2*a(s,1)));
%         lik=sum(log(normpdf(ctraj(s,2:end)',mean_h',sqrt(var_h))));
%         jac=log(1/a(s,1));
%         post=prior+lik;
%         theta=[log(a(s,1))];
%         R1=mvnrnd(theta,V1(:,:,s));
%         a_star=exp(R1(1,1));
%         prior_star=log(unifpdf(a_star,0,100));
%         mean_h_star=exp(-a_star)*ctraj(s,1:T-1);
%         var_h_star=tau(s,1)*((1-exp(-2*a_star))/(2*a_star));
%         lik_star=sum(log(normpdf(ctraj(s,2:end)',mean_h_star',sqrt(var_h_star))));
%         post_star=prior_star+lik_star;
%         jac_star=log(1/a_star);
%         r1 = exp(post_star-post+jac-jac_star);
%         C1 = min(1,r1);
%         if A1<=C1
%            a(s,1)=a_star;
%         end
%         thetasave(i,:,s)=theta;
%         if i>50
%            V1(:,:,s)=cov(thetasave(1:i,:,s));
%         end
%     end
%     
%     for s=dim_y+1:dim_y+num_fact
%         v1=v0+(T-1);
%         s1=s0+((2*a(s,1))/(1-exp(-2*a(s,1))))*sum((ctraj(s,2:end)-exp(-a(s,1))*ctraj(s,1:T-1)).^2);
%         tau(s,1)=1./random('gam',v1/2,2/s1);
%         
%         A1=rand();
%         prior=log(unifpdf(a(s,1),0,100));
%         mean_h=exp(-a(s,1))*ctraj(s,1:T-1);
%         var_h=tau(s,1)*((1-exp(-2*a(s,1)))/(2*a(s,1)));
%         lik=sum(log(normpdf(ctraj(s,2:end)',mean_h',sqrt(var_h))));
%         jac=log(1/a(s,1));
%         post=prior+lik;
%         theta=[log(a(s,1))];
%         R1=mvnrnd(theta,V2(:,:,s-dim_y));
%         a_star=exp(R1(1,1));
%         prior_star=log(unifpdf(a_star,0,100));
%         mean_h_star=exp(-a_star)*ctraj(s,1:T-1);
%         var_h_star=tau(s,1)*((1-exp(-2*a_star))/(2*a_star));
%         lik_star=sum(log(normpdf(ctraj(s,2:end)',mean_h_star',sqrt(var_h_star))));
%         post_star=prior_star+lik_star;
%         jac_star=log(1/a_star);
%         r1 = exp(post_star-post+jac-jac_star);
%         C1 = min(1,r1);
%         if A1<=C1
%            a(s,1)=a_star;
%         end
%         thetasave2(i,:,s)=theta;
%         if i>50
%            V2(:,:,s)=cov(thetasave2(1:i,:,s));
%         end
%         
%     end

% load('OU_sim.mat');
% a=a_true;
% tau=tau_true;
% burn=1000;
% nit=10000;
% iter=burn+nit;
% y=y_true;
% T=length(y);
% N=500;
% v0=10;
% s0=1;
% 
% [X,W,A,~]=smc_diffusion_OU(y,T,N,a,tau);
% ctraj=ancestor_trace(X,W,A);
% accept=0;
% 
% for i=1:iter
%     i
%     a
%     tau
%     [X,W,A,~]=csmc_diffusion_OU(y,T,N,a,tau,ctraj);
%     ctraj=ancestor_trace(X,W,A);
%     
%     A1=rand();
%     prior=log(unifpdf(a,0,100));
%     mean_h=exp(-a)*ctraj(1,1:T-1);
%     var_h=tau*((1-exp(-2*a))/(2*a));
%     lik=sum(log(normpdf(ctraj(1,2:end)',mean_h',sqrt(var_h))));
%     jac=log(1/a);
%     post=prior+lik;
%     theta=[log(a)];
%     if i<50
%        cov_mat=0.01*eye(1);
%     else
%        cov_mat=cov([log(a_store)]);
%     end
%     R1=mvnrnd(theta,cov_mat);
%     a_star=exp(R1(1,1));
%     
%     prior_star=log(unifpdf(a_star,0,100));
%     mean_h_star=exp(-a_star)*ctraj(1,1:T-1);
%     var_h_star=tau*((1-exp(-2*a_star))/(2*a_star));
%     lik_star=sum(log(normpdf(ctraj(1,2:end)',mean_h_star',sqrt(var_h_star))));
%     post_star=prior_star+lik_star;
%     jac_star=log(1/a_star);
%     r1 = exp(post_star-post+jac-jac_star);
%     C1 = min(1,r1);
%     if A1<=C1
%        a=a_star;
%     end
%     
%      v1=v0+(T-1);
%      s1=s0+((2*a)/(1-exp(-2*a)))*sum((ctraj(1,2:end)-exp(-a)*ctraj(1,1:T-1)).^2);
%      tau=1./random('gam',v1/2,2/s1);
%     
%     a_store(i,1) = a;
%     tau_store(i,1) = tau;
%     ctraj_store(i,:) = ctraj;
% end
% save('PG_OU2.mat','a_store','tau_store','ctraj_store');

%     A1=rand();
%     theta=[log(a),log(tau)];
%     if i<50
%         cov_mat=0.001*eye(2);
%     else
%         cov_mat=cov([log(a_store),log(tau_store)]);
%     end
%     R1=mvnrnd(theta,cov_mat);
%     a_star=exp(R1(1,1));
%     tau_star=exp(R1(1,2));
%     prior_star=log(unifpdf(a_star,0,100))+log_IG_PDF_used(tau_star,v0/2,s0/2);
%     [~,~,~,lik_star]=smc_diffusion_OU(y,T,N,a_star,tau_star);
%     post_star=prior_star+lik_star;
%     jac_star=log(1/a_star)+log(1/tau_star);
%     r1 = exp(post_star-post+jac-jac_star);
%     C1 = min(1,r1);
%     if A1<=C1
%        a=a_star;
%        tau=tau_star;
%        post=post_star;
%        accept=accept+1;
%        jac=jac_star;
%     end
    
% a=a_true;
% b=b_true;
% s=s_true;
% burn=1000;
% nit=10000;
% s_tot=burn+nit;
% T=length(r_true);
% y=r_true;
% %beta_store=zeros(nit,1);
% %alpha_store=zeros(nit,1);
% %sigma_store=zeros(nit,1);
% N=2500;
% M=10;
% delta=1/12;
% 
% prior=log(unifpdf(a,0,100))+log(unifpdf(s,0,100))+log(unifpdf(b,0,1));
% lik=lik_diffusion3(y,T,M,N,delta,a,b,s);
% %parfor i=1:100
% %    i
% %[lik(i,1)]=lik_diffusion3(y,T,M,N,delta,a,b,s);
% %end
% post=prior+lik;
% jac=log(1/a)+log(1/s)+log((1/b)+(1/(1-b)));
% accept=0;
% for i=1:s_tot
%     i
%     a
%     b
%     s
%     A1=rand();
%     theta=[log(a),log(s),logit_inverse(b)];
%     if i<100
%        cov_mat=0.001*eye(3);
%     else
%        %cov_mat=var(log(a_store));
%         cov_mat=cov([log(a_store),log(s_store),logit_inverse(b_store)]);
%     end
%     R1=mvnrnd(theta,cov_mat);
%     a_star=exp(R1(1,1));
%     s_star=exp(R1(1,2));
%     b_star=logit_cdf_com(R1(1,3));
%     prior_star=log(unifpdf(a_star,0,100))+log(unifpdf(s_star,0,100))+log(unifpdf(b_star,0,1));
%     lik_star=lik_diffusion3(y,T,M,N,delta,a_star,b_star,s_star);
%     post_star=prior_star+lik_star;
%     jac_star=log(1/a_star)+log(1/s_star)+log((1/b_star)+(1/(1-b_star)));
%     r1 = exp(post_star-post+jac-jac_star);
%     C1 = min(1,r1);
%     if A1<=C1
%          a=a_star;
%          s=s_star;
%          b=b_star;
%          post=post_star;
%          accept=accept+1;
%          jac=jac_star;
%     end
%     a_store(i,1) = a;
%     s_store(i,1) = s;
%     b_store(i,1) = b;
%     
% end


% load('signal_noise.mat');
% load('cov_mat.mat');
% burn=1000;
% nit=10000;
% s=burn+nit;
% T=length(y);
% sigma=0.5;
% phi=0.95;
% tau=0.5;
% tau_store=zeros(nit,1);
% sigma_store=zeros(nit,1);
% phi_store=zeros(nit,1);
% 
% N=40;
% prior=logcauchy(tau)+logcauchy(sigma)+log(unifpdf(phi,0,1));
% %for i=1:50
% u1=randn(T+1,N+1);
% lik=sirfilter_corr(tau,T,N,sigma,phi,y,u1,1);
% %end
% post=prior+lik;
% jac=log(1/tau)+log(1/sigma)+log((1/phi)+(1/(1-phi)));
% pstar=0.30;
% alpha=-norminv(pstar/2);
% n0=round(5/(pstar*(1-pstar)));
% dim_sigma1=3;
% sigma_1=1;
% sigma2_1=sigma_1^2;
% sigma_vec_1=sigma_1;
% accept=0;
% rho=0.9999;
% for i=1:s
%     i
%     tau
% 
%     A1=rand();
%     u1_star=rho*u1+sqrt(1-rho^2)*randn(T+1,N+1);
%     theta=[log(tau),log(sigma),logit_inverse(phi)];
%     R1=mvnrnd(theta,sigma2_1*cov_mat_prop1);
%     tau_star=exp(R1(1,1));
%     sigma_star=exp(R1(1,2));
%     phi_star=logit_cdf_com(R1(1,3));
%     prior_star=logcauchy(tau_star)+logcauchy(sigma_star)+log(unifpdf(phi_star,0,1));
%     lik_star=sirfilter_corr(tau_star,T,N,sigma_star,phi_star,y,u1_star,1);
%     post_star=prior_star+lik_star;
%     jac_star=log(1/tau_star)+log(1/sigma_star)+log((1/phi_star)+(1/(1-phi_star)));
%     r1 = exp(post_star-post+jac-jac_star);
%     
%     C1 = min(1,r1);
%     
%     if A1<=C1
%         tau=tau_star;
%         sigma=sigma_star;
%         phi=phi_star;
%         u1=u1_star;
%         post=post_star;
%         accept=accept+1;
%         jac=jac_star;
%     else
%         tau=tau;
%         sigma=sigma;
%         phi=phi;
%         post=post;
%         jac=jac;
%         u1=u1;
%     end
%     
%      if i>burn
%          tau_store(i-burn,1) = tau;
%          sigma_store(i-burn,1)= sigma;
%          phi_store(i-burn,1) = phi;
%     end
%          tau_tot(i,1) = tau;
%          sigma_tot(i,1) = sigma;
%          phi_tot(i,1) = phi;
%     if i>n0
%          sigma_1=update_sigma(sigma2_1,C1,pstar,i,dim_sigma1);
%          sigma2_1=sigma_1^2;
%          sigma_vec_1=[sigma_vec_1;sigma_1];
%      end
% end

%     if i<=50
%         diag_cov=[0.05;0.05];
%         cov_mat_prop1=diag(diag_cov);
%     else
%         theta_sub_tot=[log(tau_tot),log(sigma_tot)];
%         cov_mat_prop1=cov(theta_sub_tot);
%     end


