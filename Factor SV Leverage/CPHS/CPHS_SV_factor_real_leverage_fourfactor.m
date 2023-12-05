%PMMH pf
%load('DATA_REAL30_SV_general.mat');
%parpool(28)
%matlabpool open local 12

load('dataindustrystockreal3000.mat');
dim_y=26; %the dimensions of stock returns.
num_fact=4; %the number of factors
%initial values of the parameters
phi=0.98*ones(dim_y+num_fact,1);
tau=0.1*ones(dim_y+num_fact,1);
mu=0.01*ones(dim_y,1);
rho=-0.1*ones(dim_y,1);
%B(:,1)=rand(dim_y,1);
B(:,1)=rand(dim_y,1);
B(:,2)=[0;rand(dim_y-1,1)];
B(:,3)=[0;0;rand(dim_y-2,1)];
B(:,4)=[0;0;0;rand(dim_y-3,1)];
%the number of MCMC iterations
burn=1000;
nit=11000;
iter=burn+nit;
N=100;%the number of particles
T=3001;%the length of time series
y=y_order_demean_used;
%initial scale and covariance matrix of the random walk proposals.
D1=2;
for i=1:dim_y
    scale_idio(i,1)=1;
    V1(:,:,i)=0.01*eye(D1);
end
target_accept=0.20;
D2=1;
for i=1:num_fact
    scale_factor(i,1)=1;
    V2(:,:,i)=0.01;
end
%initial values of the latent factors
fact_score=randn(num_fact,T);
%initial values of the log-volatilities
parfor s=1:dim_y
  [X,W,A,~]=smc_mu_factorSV_leverage(y(s,:),B(s,:),fact_score,phi(s,1),tau(s,1),mu(s,1),rho(s,1),N,T,s,dim_y);
  ctraj(s,:)=ancestor_trace(X,W,A);
end
parfor s=dim_y+1:dim_y+num_fact
  [X,W,A,~]=smc_mu_factorSV(0,0,fact_score,phi(s,1),tau(s,1),0,N,T,s,dim_y);
  ctraj(s,:)=ancestor_trace(X,W,A);
end
%prior hyperparameters
B0_deep=10^12;
hp_sig2=1;
v0=10;
s0=0.25;

a_phi=100;
b_phi=1.5;

accept1=zeros(num_fact,1);
accept2=zeros(dim_y+num_fact,1);
hp_mu=4;

ctraj1_store=zeros(1,T);
ctraj2_store=zeros(1,T);
ctraj3_store=zeros(1,T);
ctraj4_store=zeros(1,T);
ctraj5_store=zeros(1,T);
ctraj6_store=zeros(1,T);
ctraj7_store=zeros(1,T);
ctraj8_store=zeros(1,T);
ctraj9_store=zeros(1,T);
ctraj10_store=zeros(1,T);
ctraj11_store=zeros(1,T);
ctraj12_store=zeros(1,T);
ctraj13_store=zeros(1,T);
ctraj14_store=zeros(1,T);
ctraj15_store=zeros(1,T);
ctraj16_store=zeros(1,T);
ctraj17_store=zeros(1,T);
ctraj18_store=zeros(1,T);
ctraj19_store=zeros(1,T);
ctraj20_store=zeros(1,T);
ctraj21_store=zeros(1,T);
ctraj22_store=zeros(1,T);
ctraj23_store=zeros(1,T);
ctraj24_store=zeros(1,T);
ctraj25_store=zeros(1,T);
ctraj26_store=zeros(1,T);

ctrajf1_store=zeros(1,T);
ctrajf2_store=zeros(1,T);
ctrajf3_store=zeros(1,T);
ctrajf4_store=zeros(1,T);
ctrajf5_store=zeros(1,T);

for i=1:iter
    
        
	%sampling factor loadings in Gibbs steps 
    for s=1:dim_y
        s_bar=min(s,num_fact);
        F=fact_score(1:s_bar,:)';
        arg_Vi_inv=1./(exp(ctraj(s,:)));
        Vi_inv=diag(arg_Vi_inv);
        B_var=inv((F'*Vi_inv*F)+eye(s_bar));
        [B_var]=jitChol(B_var);
        chol_B=chol(B_var);
        B_mean=B_var*F'*((Vi_inv*y(s,:)'));
        B(s,1:s_bar)=(mvnrnd(B_mean',chol_B'*chol_B))';        
    end
        
    
     %deep interweaving algorithm of Kastner et al (2017)
     k1=dim_y;
     for s=1:num_fact
       ind=find(abs(B(:,s))==max(abs(B(:,s))));
       max_B=(B(s,s));
       B_star=B(s:dim_y,s)./max_B;
       fact_star=max_B.*fact_score(s,:);
       ctraj_trans=ctraj(dim_y+s,:)+2*log(abs(max_B));
       mu_fact_old=log(max_B^2);
       mean_prop=(sum(ctraj_trans(1,2:T-1))+((ctraj_trans(1,T)-phi(dim_y+s,1)*ctraj_trans(1,1))./(1-phi(dim_y+s,1))))./(T-1+1/B0_deep);
       var_prop=(tau(dim_y+s,1)/((1-phi(dim_y+s,1))^2))/(T-1+1/B0_deep);
       mu_fact_new=normrnd(mean_prop,sqrt(var_prop));
       A2=rand();
       comp1_old=logmvnpdf(B_star',zeros(1,k1),hp_sig2*exp(-mu_fact_old)*eye(k1));
       comp2_old=logmvnpdf(ctraj_trans(1,1),mu_fact_old,tau(dim_y+s,1)/(1-phi(dim_y+s,1)^2));       
       comp3_old=log(exp(mu_fact_old/2-exp(mu_fact_old)/(2*hp_sig2)));
       comp4_old=logmvnpdf(mu_fact_old,0,B0_deep*tau(dim_y+s,1)/((1-phi(dim_y+s,1))^2));
       comp1_new=logmvnpdf(B_star',zeros(1,k1),hp_sig2*exp(-mu_fact_new)*eye(k1));
       comp2_new=logmvnpdf(ctraj_trans(1,1),mu_fact_new,tau(dim_y+s,1)/(1-phi(dim_y+s,1)^2));
       comp3_new=log(exp(mu_fact_new/2-exp(mu_fact_new)/(2*hp_sig2)));
       comp4_new=logmvnpdf(mu_fact_new,0,B0_deep*tau(dim_y+s,1)/((1-phi(dim_y+s,1))^2));
       R2=exp(comp1_new+comp2_new+comp3_new-comp1_old-comp2_old-comp3_old+comp4_old-comp4_new);
       C2=min(1,R2);
       if A2<=C2
          lam1=exp(mu_fact_new/2);
          accept1(s,1)=accept1(s,1)+1;   
       else
          lam1=max_B;
       end     
       B(:,s)=(lam1/max_B).*B(:,s);
       fact_score(s,:)=(max_B/lam1).*fact_score(s,:);
       ctraj(dim_y+s,:)=ctraj(dim_y+s,:)+2*log(abs(max_B/lam1));
       k1=k1-1;
     end

     %sampling the latent factor ft, for t=1,..,T.
    for t=1:T
        ctraj_idio_temp(:,t)=ctraj(1:dim_y,t);
        ctraj_fact_temp(:,t)=ctraj(dim_y+1:dim_y+num_fact,t);
        Vt_inv=diag(1./(exp(ctraj_idio_temp(:,t))));
        Dt_inv=diag(exp(-ctraj_fact_temp(:,t)));
        var_ft=inv((B'*Vt_inv*B)+Dt_inv);
        [var_ft]=jitChol(var_ft);
        chol_var_ft=chol(var_ft);
        mean_ft=var_ft*B'*(Vt_inv*y(:,t));
        fact_score(:,t)=mvnrnd(mean_ft,chol_var_ft'*chol_var_ft);      
    end
    %sampling the SV parameters mu and phi in Gibbs steps
    for s=1:dim_y
         eps=(y(s,1:T-1)-B(s,:)*fact_score(:,1:T-1)).*exp(-ctraj(s,1:T-1)./2);         

         mu_var=(tau(s,1)*(1-rho(s,1)^2))/((1-rho(s,1)^2)*(1-phi(s,1)^2)+(T-1)*(1-phi(s,1))^2);
         mu_mean=(mu_var/(tau(s,1)*(1-rho(s,1)^2)))*((1-rho(s,1)^2)*(1-phi(s,1)^2)*(ctraj(s,1))+(1-phi(s,1)).*sum(ctraj(s,2:T)-phi(s,1).*(ctraj(s,1:T-1))-rho(s,1).*sqrt(tau(s,1)).*eps));
         mu(s,1)=normrnd(mu_mean,sqrt(mu_var));
          
         %phi_var=(tau(s,1)*(1-rho(s,1)^2))/(rho(s,1)^2*(ctraj(s,1)-mu(s,1)).^2+sum(ctraj(s,2:T-1)-mu(s,1)).^2);
         %phi_mean=sum(((ctraj(s,2:T)-mu(s,1))-sqrt(tau(s,1)).*rho(s,1).*eps).*(ctraj(s,1:T-1)-mu(s,1)))./(rho(s,1)^2*(ctraj(s,1)-mu(s,1)).^2+sum(ctraj(s,2:T-1)-mu(s,1)).^2);
         
         phi_var=(tau(s,1)*(1-rho(s,1)^2))/(sum((ctraj(s,1:T-1)-mu(s,1)).^2)-((ctraj(s,1)-mu(s,1))^2)*(1-rho(s,1)^2));
         phi_mean=sum((ctraj(s,2:T)-mu(s,1)).*(ctraj(s,1:T-1)-mu(s,1))-rho(s,1).*sqrt(tau(s,1)).*(ctraj(s,1:T-1)-mu(s,1)).*eps)./(sum((ctraj(s,1:T-1)-mu(s,1)).^2)-((ctraj(s,1)-mu(s,1))^2)*(1-rho(s,1)^2));
         if phi_mean>0.999
         phi_mean=0.999;
         end
         phi_star=normt_rnd(phi_mean,phi_var,0,0.999);
         prior=log(betapdf((1+phi(s,1))/2,a_phi,b_phi));
         prior_star=log(betapdf((1+phi_star)/2,a_phi,b_phi));
         num_phi=log(sqrt(1-phi_star^2));
         den_phi=log(sqrt(1-phi(s,1)^2));
         r2=exp(num_phi-den_phi+prior_star-prior);
         C2=min(1,r2);
         A2=rand();
         if A2<=C2
            phi(s,1)=phi_star;
         end
    end
    
    for s=dim_y+1:dim_y+num_fact
       phi_var=1./(sum((ctraj(s,2:T-1).^2)./tau(s,1)));
       phi_mean=phi_var.*(sum((ctraj(s,2:T).*ctraj(s,1:T-1))./tau(s,1)));
       if phi_mean>1
          phi_mean=1;
       end
       phi_star=normt_rnd(phi_mean,phi_var,0,0.99999);
       num_phi=log(sqrt(1-phi_star^2));
       den_phi=log(sqrt(1-phi(s,1)^2));
       prior=log(betapdf((1+phi(s,1))/2,a_phi,b_phi));
       prior_star=log(betapdf((1+phi_star)/2,a_phi,b_phi));
       r2=exp(num_phi-den_phi+prior_star-prior);
       C2=min(1,r2);
       A2=rand();
       if A2<=C2
          phi(s,1)=phi_star;
       end
    end
    %sampling tau2 and rho in pseudo marginal steps
    parfor s=1:dim_y
        [u1]=obtain_random_numbers_leverage(phi(s,1),tau(s,1),mu(s,1),rho(s,1),ctraj(s,:),N,T,y(s,:),B(s,:),fact_score); 
        [X,W,A,lik,u1_res_rand]=csmc_mu_factorSV_corr_leverage(y(s,:),B(s,:),fact_score,phi(s,1),tau(s,1),mu(s,1),rho(s,1),N,T,ctraj(s,:),u1,s,dim_y);        
        prior=logcauchy(tau(s,1));
        post=prior+lik;
        jac=log(1/tau(s,1));
        theta=[log(tau(s,1)),atanh(rho(s,1))];
        R1=mvnrnd(theta,scale_idio(s,1).*V1(:,:,s));
        tau_star=exp(R1(1,1));
        rho_star=tanh(R1(1,2));
        [X_star,W_star,A_star,lik_star]=smc_mu_factorSV_corr_leverage(y(s,:),B(s,:),fact_score,phi(s,1),tau_star,mu(s,1),rho_star,N,T,u1,u1_res_rand,s,dim_y);
        prior_star=logcauchy(tau_star);
        %prior_star=log_IG_PDF_used(tau_star,v0/2,s0/2);
        post_star=prior_star+lik_star;
        jac_star=log(1/tau_star);
        r1=exp(post_star-post+jac-jac_star);
        C1 = min(1,r1);
        A1=rand();
        if A1<=C1
           tau(s,1)=tau_star; 
           rho(s,1)=rho_star;
           X=X_star;
           W=W_star;
           A=A_star;
           accept2(s,1)=accept2(s,1)+1; 
        end
        thetasave{s,1}(i,:)=[log(tau(s,1)),atanh(rho(s,1))];
        ctraj(s,:)=backward_simulation_SV_leverage(X,W,A,phi(s,1),tau(s,1),mu(s,1),rho(s,1),y(s,:),B(s,:),fact_score);
        if i>50
           scale_idio(s,1)=update_sigma(scale_idio(s,1),C1,target_accept,i,2); 
        end
    end
    
    for s=1:dim_y
         if i>50
            V1(:,:,s)=cov(thetasave{s,1}(1:i,:));
            V1(:,:,s)=jitChol(V1(:,:,s));
         end
    end
 
    parfor s=dim_y+1:dim_y+num_fact
       [u1]=obtain_random_numbers(phi(s,1),tau(s,1),0,ctraj(s,:),N,T);
       [X,W,A,lik,u1_res_rand]=csmc_mu_factorSV_corr(0,0,fact_score,phi(s,1),tau(s,1),0,N,T,ctraj(s,:),u1,s,dim_y);
       %prior=log_IG_PDF_used(tau(s,1),v0/2,s0/2);
       prior=logcauchy(tau(s,1));
       post=prior+lik;
       jac=log(1/tau(s,1));
       theta=log(tau(s,1));
       R1=mvnrnd(theta,1.5.*V2(:,:,s-dim_y));
       tau_star=exp(R1(1,1));
       [X_star,W_star,A_star,lik_star]=smc_mu_factorSV_corr(0,0,fact_score,phi(s,1),tau_star,0,N,T,u1,u1_res_rand,s,dim_y);
       prior_star=logcauchy(tau_star);
       %prior_star=log_IG_PDF_used(tau_star,v0/2,s0/2);
       post_star=prior_star+lik_star;
       jac_star=log(1/tau_star);
       r1 = exp(post_star-post+jac-jac_star);
       C1 = min(1,r1);
       A1=rand();
       if A1<=C1
          tau(s,1)=tau_star; 
          X=X_star;
          W=W_star;
          A=A_star;
          accept2(s,1)=accept2(s,1)+1; 
       end
       thetasave2(i,s-dim_y)=log(tau(s,1));
       ctraj(s,:)=backward_simulation_SV(X,W,A,phi(s,1),tau(s,1),0); 
       if i>50
           scale_factor(s-dim_y,1)=update_sigma(scale_factor(s-dim_y,1),C1,target_accept,i,1); 
        end
    end
    
    for s=dim_y+1:dim_y+num_fact
         if i>50
            V2(:,:,s-dim_y)=cov(thetasave2(1:i,s-dim_y));
            V2(:,:,s-dim_y)=jitChol(V2(:,:,s-dim_y)); 
         end 
    end
    
    %computing DIC
    %computing first term and second term
    for t=1:T
        
         %computing first term
         arg_Vt=exp(ctraj(1:dim_y,t));
         Vt=diag(arg_Vt);
         Dt_arg=exp(ctraj(dim_y+1:dim_y+num_fact,t));
         Dt=diag(Dt_arg);
         var_arg=B*Dt*B'+Vt;
         chol_var_arg=chol(var_arg);
         llh_ind_DIC(t,1)=logmvnpdf(y(:,t)',zeros(1,dim_y),chol_var_arg'*chol_var_arg); 
     
         %computing second term
         if t==1
            llh_ctraj_idio_ind=log(normpdf(ctraj(1:dim_y,1),mu(1:dim_y),sqrt(tau(1:dim_y)./(1-phi(1:dim_y).^2))));
            llh_ctraj_factor_ind=log(normpdf(ctraj(dim_y+1:dim_y+num_fact,1),0,sqrt(tau(dim_y+1:dim_y+num_fact)./(1-phi(dim_y+1:dim_y+num_fact).^2))));
         else
            epsilon_tmin1=(y(1:dim_y,t-1)-B*fact_score(:,t-1)).*exp(-ctraj(1:dim_y,1)./2); 
            mean_com_idio=mu(1:dim_y,1)+phi(1:dim_y,1).*(ctraj(1:dim_y,t-1)-mu(1:dim_y,1))+rho(1:dim_y,1).*sqrt(tau(1:dim_y,1)).*epsilon_tmin1;
            var_com_idio=(1-rho(1:dim_y,1).^2).*tau(1:dim_y,1);
            llh_ctraj_idio_ind=llh_ctraj_idio_ind+log(normpdf(ctraj(1:dim_y,t),mean_com_idio,sqrt(var_com_idio)));
            mean_com_factor=phi(dim_y+1:dim_y+num_fact,1).*ctraj(dim_y+1:dim_y+num_fact,t-1);
            var_com_factor=tau(dim_y+1:dim_y+num_fact,1);
            llh_ctraj_factor_ind=llh_ctraj_factor_ind+log(normpdf(ctraj(dim_y+1:dim_y+num_fact,t),mean_com_factor,sqrt(var_com_factor)));
         end 
    end
    %computing prior
     col_B=[B(:,1);B(2:end,2);B(3:end,3)];
     prior_B=sum(log(normpdf(col_B,0,1)));
     prior_tau=sum(logcauchy(tau));
     %prior_tau=sum(log_IG_PDF_used(tau,v0/2,s0/2));
     prior_phi=sum(log(betapdf((1+phi)./2,a_phi,b_phi)));
     prior_tot=prior_B+prior_tau+prior_phi;
    
     llh_ctraj_idio=sum(llh_ctraj_idio_ind);
     llh_ctraj_factor=sum(llh_ctraj_factor_ind);
     
     llh_DIC=sum(llh_ind_DIC);   
     llh_store(i,1)=llh_DIC;
     llh_post(i,1)=llh_DIC+llh_ctraj_idio+llh_ctraj_factor+prior_tot;
    
     if i>burn

        ctraj1_store=ctraj1_store+ctraj(1,:);    
        ctraj2_store=ctraj2_store+ctraj(2,:);    
        ctraj3_store=ctraj3_store+ctraj(3,:);    
        ctraj4_store=ctraj4_store+ctraj(4,:);    
        ctraj5_store=ctraj5_store+ctraj(5,:);    
        ctraj6_store=ctraj6_store+ctraj(6,:);    
        ctraj7_store=ctraj7_store+ctraj(7,:);    
        ctraj8_store=ctraj8_store+ctraj(8,:);    
        ctraj9_store=ctraj9_store+ctraj(9,:);    
        ctraj10_store=ctraj10_store+ctraj(10,:);    
        ctraj11_store=ctraj11_store+ctraj(11,:);    
        ctraj12_store=ctraj12_store+ctraj(12,:);    
        ctraj13_store=ctraj13_store+ctraj(13,:);    
        ctraj14_store=ctraj14_store+ctraj(14,:);    
        ctraj15_store=ctraj15_store+ctraj(15,:);    
        ctraj16_store=ctraj16_store+ctraj(16,:);    
        ctraj17_store=ctraj17_store+ctraj(17,:);    
        ctraj18_store=ctraj18_store+ctraj(18,:);    
        ctraj19_store=ctraj19_store+ctraj(19,:);    
        ctraj20_store=ctraj20_store+ctraj(20,:);    
        ctraj21_store=ctraj21_store+ctraj(21,:);    
        ctraj22_store=ctraj22_store+ctraj(22,:);    
        ctraj23_store=ctraj23_store+ctraj(23,:);    
        ctraj24_store=ctraj24_store+ctraj(24,:);    
        ctraj25_store=ctraj25_store+ctraj(25,:);    
        ctraj26_store=ctraj26_store+ctraj(26,:);         
        ctrajf1_store=ctrajf1_store+ctraj(27,:);    
        ctrajf2_store=ctrajf2_store+ctraj(28,:);    
        ctrajf3_store=ctrajf3_store+ctraj(29,:);    
        ctrajf4_store=ctrajf4_store+ctraj(30,:);     
   
    end
     
     

    Post.B1(i,:)=B(:,1)';
    Post.B2(i,:)=B(:,2)';
    Post.B3(i,:)=B(:,3)';
    Post.B4(i,:)=B(:,4)';
    Post.tau(i,:)=tau';
    Post.phi(i,:)=phi';
    Post.mu(i,:)=mu';
    Post.rho(i,:)=rho';
    id=(1:1:100)*20;
    Post.ctraj1(i,:)=ctraj(1,id);
    Post.ctraj2(i,:)=ctraj(2,id);
    Post.ctraj3(i,:)=ctraj(3,id);
    Post.ctraj4(i,:)=ctraj(4,id);
    Post.ctraj5(i,:)=ctraj(5,id);
    Post.ctraj6(i,:)=ctraj(6,id);
    Post.ctraj7(i,:)=ctraj(7,id);
    Post.ctraj8(i,:)=ctraj(8,id);
    Post.ctraj9(i,:)=ctraj(9,id);
    Post.ctraj10(i,:)=ctraj(10,id);
    Post.ctraj11(i,:)=ctraj(11,id);
    Post.ctraj12(i,:)=ctraj(12,id);
    Post.ctraj13(i,:)=ctraj(13,id);
    Post.ctraj14(i,:)=ctraj(14,id);
    Post.ctraj15(i,:)=ctraj(15,id);
    Post.ctraj16(i,:)=ctraj(16,id);
    Post.ctraj17(i,:)=ctraj(17,id);
    Post.ctraj18(i,:)=ctraj(18,id);
    Post.ctraj19(i,:)=ctraj(19,id);
    Post.ctraj20(i,:)=ctraj(20,id);
    Post.ctraj21(i,:)=ctraj(21,id);
    Post.ctraj22(i,:)=ctraj(22,id);
    Post.ctraj23(i,:)=ctraj(23,id);
    Post.ctraj24(i,:)=ctraj(24,id);
    Post.ctraj25(i,:)=ctraj(25,id);
    Post.ctraj26(i,:)=ctraj(26,id);
    Post.ctrajfactor1(i,:)=ctraj(27,id);
    Post.ctrajfactor2(i,:)=ctraj(28,id);
    Post.ctrajfactor3(i,:)=ctraj(29,id);
    Post.ctrajfactor4(i,:)=ctraj(30,id);
    
    Post.ctraj1sum=ctraj1_store;
    Post.ctraj2sum=ctraj2_store;
    Post.ctraj3sum=ctraj3_store;
    Post.ctraj4sum=ctraj4_store;
    Post.ctraj5sum=ctraj5_store;
    Post.ctraj6sum=ctraj6_store;
    Post.ctraj7sum=ctraj7_store;
    Post.ctraj8sum=ctraj8_store;
    Post.ctraj9sum=ctraj9_store;
    Post.ctraj10sum=ctraj10_store;
    Post.ctraj11sum=ctraj11_store;
    Post.ctraj12sum=ctraj12_store;
    Post.ctraj13sum=ctraj13_store;
    Post.ctraj14sum=ctraj14_store;
    Post.ctraj15sum=ctraj15_store;
    Post.ctraj16sum=ctraj16_store;
    Post.ctraj17sum=ctraj17_store;
    Post.ctraj18sum=ctraj18_store;
    Post.ctraj19sum=ctraj19_store;
    Post.ctraj20sum=ctraj20_store;
    Post.ctraj21sum=ctraj21_store;
    Post.ctraj22sum=ctraj22_store;
    Post.ctraj23sum=ctraj23_store;
    Post.ctraj24sum=ctraj24_store;
    Post.ctraj25sum=ctraj25_store;
    Post.ctraj26sum=ctraj26_store;
    Post.ctrajf1sum=ctrajf1_store;
    Post.ctrajf2sum=ctrajf2_store;
    Post.ctrajf3sum=ctrajf3_store;
    Post.ctrajf4sum=ctrajf4_store;
    
     if i==250 | i==500 | i==1000 | i==2000 | i==3000 | i==4000 | i==5000 | i==6000 | i==7000 | i==10000 | i==12000 | i==14000
        save('/short/jz21/dg2271/corr_mix/factor_SV/corr_factor_real_leverage4factorN20.mat','Post','llh_store','llh_post');
     end
     
end
save('/short/jz21/dg2271/corr_mix/factor_SV/corr_factor_real_leverage4factorN20.mat','Post','llh_store','llh_post');

