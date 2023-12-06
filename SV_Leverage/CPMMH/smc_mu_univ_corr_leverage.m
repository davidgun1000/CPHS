function [particles,w,indx,sir_llh]=smc_mu_univ_corr_leverage(y,phi,tau,mu,rho,N,T,u1,u1_res)
   t=1; 
   particles=zeros(1,N,T); 
   w=zeros(T,N);
   indx=zeros(T,N); 
   particles(:,1:N,t)=sqrt(tau/(1-phi^2))*u1(:,t)'+mu;
   logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:));
   for t=2:T
       [indx(t,:)]=rs_multinomial_corr2_PMMH(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
       eps_t=(y(1,t-1)).*exp(-particles(:,indx(t,1:N),t-1)./2);
       particles(:,1:N,t)=mu+phi*(particles(:,indx(t,1:N),t-1)-mu)+rho*sqrt(tau).*eps_t+sqrt(1-rho^2).*sqrt(tau).*u1(:,t)';
       logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
       w(t,:)=exp(logw-max(logw));
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end

end







