function [particles,w,indx,sir_llh]=csmc_mu_SV_Fearnhead_leverage(y,pseudo_obs,pseudo_hyper,prior_hyper,current_param_trans,N,T,ctraj)

   t=1;
   particles=zeros(5,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(:,1,:)=ctraj;
   current_param.phi=logit_cdf_com(current_param_trans.phi);
   current_param.mu=current_param_trans.mu;
   current_param.tau=exp(current_param_trans.tau);
   current_param.rho=tanh(current_param_trans.rho);
   
   prior_phi_trans_current=log(betapdf((current_param.phi+1)/2,prior_hyper.a0,prior_hyper.b0))+...
        log(exp(current_param_trans.phi))-log((1+exp(current_param_trans.phi)).^2);
     phi_trans_prop=normrnd(pseudo_obs.phi_trans,sqrt(pseudo_hyper.phi_trans),1,N-1);
     prop_param.phi=logit_cdf_com(phi_trans_prop);
     prior_phi_trans_prop=log(betapdf((prop_param.phi+1)/2,prior_hyper.a0,prior_hyper.b0))+...
        log(exp(phi_trans_prop))-log((1+exp(phi_trans_prop)).^2);
     r2=exp(prior_phi_trans_prop-prior_phi_trans_current);
     C2=min(1,r2);
     A2=rand(1,N-1);
     id=(A2<=C2);
     temp_phi(1,id)=prop_param.phi(1,id);
     id=1-id;
     id=logical(id);
     temp_phi(1,id)=current_param.phi;
     particles(2,2:N,t)=temp_phi;
     
     particles(3,2:N,t)=normrnd(pseudo_obs.mu_trans,sqrt(pseudo_hyper.mu_trans),1,N-1);
     
     prior_tau_trans_current=-log(pi*(1+current_param.tau))+0.5*log(current_param.tau);
     tau_trans_prop=normrnd(pseudo_obs.tau_trans,sqrt(pseudo_hyper.tau_trans),1,N-1);
     prop_param.tau=exp(tau_trans_prop);
     prior_tau_trans_prop=-log(pi.*(1+prop_param.tau))+0.5.*log(prop_param.tau);
     r2=exp(prior_tau_trans_prop-prior_tau_trans_current);
     C2=min(1,r2);
     A2=rand(1,N-1);
     id=(A2<=C2);
     temp_tau(1,id)=prop_param.tau(1,id);
     id=1-id;
     id=logical(id);
     temp_tau(1,id)=current_param.tau;
     particles(4,2:N,t)=temp_tau;
     
     temp_rho=normrnd(pseudo_obs.rho_trans,sqrt(pseudo_hyper.rho_trans),1,N-1);
     particles(5,2:N,t)=tanh(temp_rho);
     
     particles(1,2:N,t)=sqrt(particles(4,2:N,t)./(1-particles(2,2:N,t).^2)).*randn(1,N-1)+particles(3,2:N,t);
     logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
     w(t,:)=exp(logw-max(logw));
     sir_llh=log(mean(w(t,:)))+max(logw);
     w(t,:)=w(t,:)./sum(w(t,:)); 
    for t=2:T
        
         [indx(t,:)]=rs_multinomial_cond(w(t-1,:));
         eps_t=(y(1,t-1)).*exp(-particles(1,indx(t,2:N),t-1)./2);
         particles(2,2:N,t)=particles(2,indx(t,2:N),t-1);
         particles(3,2:N,t)=particles(3,indx(t,2:N),t-1);
         particles(4,2:N,t)=particles(4,indx(t,2:N),t-1);
         particles(5,2:N,t)=particles(5,indx(t,2:N),t-1);
         particles(1,2:N,t)=particles(3,2:N,t)+particles(2,2:N,t).*(particles(1,indx(t,2:N),t-1)-particles(3,2:N,t))+...
             particles(5,2:N,t).*sqrt(particles(4,2:N,t)).*eps_t+sqrt(1-particles(5,2:N,t).^2).*sqrt(particles(4,2:N,t)).*randn(1,N-1);
         logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
         w(t,:)=exp(logw-max(logw)); 
         sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
         w(t,:)=w(t,:)./sum(w(t,:));
    end
 end