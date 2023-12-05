function [particles,w,indx,sir_llh]=csmc_mu_factorSV_leverage(y,B,ft,phi,tau,mu,rho,N,T,ctraj,num,dim)

if num<=dim
   t=1;
   particles=zeros(1,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(1,1,:)=ctraj;
   particles(:,2:N,1)=sqrt(tau/(1-phi^2))*randn(1,N-1)+mu;
   logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-B(1,:)*ft(:,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:)); 
   for t=2:T
       [indx(t,:)]=rs_multinomial_cond(w(t-1,:));
       eps_t=(y(1,t-1)-B(1,:)*ft(:,t-1)).*exp(-particles(:,indx(t,2:N),t-1)./2);
       particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+rho*sqrt(tau).*eps_t+sqrt(1-rho^2).*sqrt(tau).*randn(1,N-1);
       logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-B(1,:)*ft(:,t)).^2);
       w(t,:)=exp(logw-max(logw)); 
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end
 
else
    t=1;
    particles=zeros(1,N,T);
    w=zeros(T,N);
    indx=zeros(T,N);
    particles(1,1,:)=ctraj;
    particles(:,2:N,1)=sqrt(tau/(1-phi^2))*randn(1,N-1)+mu;
    logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*exp(-particles(1,:,t)).*(ft(num-dim,t).^2);
    w(t,:)=exp(logw-max(logw));
    sir_llh=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 
    for t=2:T
        [indx(t,:)]=rs_multinomial_cond(w(t-1,:));
        particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+sqrt(tau)*randn(1,N-1);
        logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*exp(-particles(1,:,t)).*(ft(num-dim,t).^2);
        w(t,:)=exp(logw-max(logw)); 
        sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
        w(t,:)=w(t,:)./sum(w(t,:));
    end
end
end