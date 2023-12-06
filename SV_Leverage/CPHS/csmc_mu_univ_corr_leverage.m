function [particles,w,indx,sir_llh,u_res]=csmc_mu_univ_corr_leverage(y,phi,tau,mu,rho,N,T,ctraj,u1)
%constrained conditional sequential monte carlo algorithm
   t=1;
   particles=zeros(1,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(1,1,:)=ctraj;
   particles(:,2:N,1)=sqrt(tau/(1-phi^2))*u1(2:end,t)'+mu;
   logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:)); 
   for t=2:T
       u=sort(rand(1,N-1));            
       [indx(t,:),u_res(:,t)]=rs_multinomial_cond_corr2(particles(:,1:N,t-1),w(t-1,:),u);
       eps_t=(y(1,t-1)).*exp(-particles(:,indx(t,2:N),t-1)./2);
       particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+rho*sqrt(tau).*eps_t+sqrt(1-rho^2).*sqrt(tau).*u1(2:end,t)';
       logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
       w(t,:)=exp(logw-max(logw)); 
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end
 end