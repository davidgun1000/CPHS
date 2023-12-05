function [particles,w,indx,sir_llh]=csmc_diffusion_GARCH_factor_SV_mean_approx(y,T,M,N,delta,a,tau,mu,B,ft,ctraj,num,dim,phi_factor,tau_factor)

if num<=dim
   
   t=1;
   h=delta/M;
   particles=zeros(1,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(1,1,:)=ctraj;
   particles(:,2:N,t)=normrnd(0,1,1,N-1);
   
   logw=-0.5*log(2*pi)-0.5.*log((exp(particles(1,:,t))))-0.5.*(1./(exp(particles(1,:,t)))).*((y(1,t)-B(1,:)*ft(:,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:));
   
   for t=2:T
       indx(t,:)=rs_systematic_cond(w(t-1,:));
       particles_res=particles(:,indx(t,2:N),t-1);
       z=particles_res;
       for i=1:M
           z=z+h.*(a.*(mu-exp(z)).*exp(-z)-(tau/2))+sqrt(tau)*sqrt(h)*randn(1,N-1);
           %z=z+h*(a*(mu-z))+sqrt(tau)*sqrt(h)*randn(1,N-1);
           z=real(z);
       end
       particles(:,2:N,t)=z;
       logw=-0.5*log(2*pi)-0.5.*log((exp(particles(1,:,t))))-0.5.*(1./(exp(particles(1,:,t)))).*((y(1,t)-B(1,:)*ft(:,t)).^2);
       w(t,:)=exp(logw-max(logw));
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)/sum(w(t,:)); 
   end

   
else

   t=1; 
   particles=zeros(1,N,T); 
   w=zeros(T,N);
   indx=zeros(T,N); 
   particles(1,1,:)=ctraj;
   particles(:,2:N,t)=normrnd(0,sqrt(tau_factor/(1-phi_factor^2)),1,N-1);
   logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((ft(num-dim,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:));
   for t=2:T
       indx(t,:)=rs_systematic_cond(w(t-1,:));
       particles(:,2:N,t)=phi_factor*particles(:,indx(t,2:N),t-1)+sqrt(tau_factor)*randn(1,N-1);
       logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((ft(num-dim,t)).^2);
       w(t,:)=exp(logw-max(logw));
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end
   
end
    



end

% t=1;
% particles=zeros(1,N,T);
% w=zeros(T,N);
% indx=zeros(T,N);
% particles(1,1,:)=ctraj;
% particles(:,2:N,t)=normrnd(0,1,1,N-1);
% logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
% w(t,:)=exp(logw-max(logw));
% sir_llh=log(mean(w(t,:)))+max(logw);
% w(t,:)=w(t,:)./sum(w(t,:));
% 
% for t=2:T
%     indx(t,:)=rs_systematic_cond(w(t-1,:));
%     particles(:,2:N,t)=exp(-a)*particles(:,indx(t,2:N),t-1)+sqrt(tau)*sqrt((1-exp(-2*a))/(2*a))*randn(1,N-1);
%     logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
%     w(t,:)=exp(logw-max(logw)); 
%     sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
%     w(t,:)=w(t,:)/sum(w(t,:)); 
% end






% tau=exp(theta_tau);
% t=1;
% particles=zeros(1,N,T);
% w=zeros(T,N);
% indx=zeros(T,N);
% particles(:,1:N,1)=sqrt(tau/(1-phi^2))*randn(1,N)+mu;
% %mean_arg=(rho/sqrt(tau)).*exp(particles(1,:,t)./2).*(particles(1,:,t)-mu+phi*mu);
% %var_arg=exp(particles(1,:,t)).*(1-(rho^2));
% %logw=-0.5*log(2*pi)-0.5.*log(var_arg)-0.5.*(1./var_arg).*((y(1,t)-mean_arg).^2);
% logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
% w(1,:)=exp(logw-max(logw));
% sir_llh=log(mean(w(1,:)))+max(logw(:,1:end));
% w(1,:)=w(1,:)./sum(w(1,:));
% 
% for t=2:T
%     indx(t,:)=rs_systematic(w(t-1,:),false);
%     particles(:,1:N,t)=mu+phi*(particles(:,indx(t,1:N),t-1)-mu)+sqrt(tau)*randn(1,N);
%     %mean_arg=(rho/sqrt(tau)).*exp(particles(1,:,t)./2).*(particles(1,:,t)-(mu+phi.*(particles(1,indx(t,1:N),t-1)-mu)));
%     %var_arg=exp(particles(1,:,t)).*(1-(rho^2));
%     %logw=-0.5*log(2*pi)-0.5.*log(var_arg)-0.5.*(1./var_arg).*((y(1,t)-mean_arg).^2);
%     logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
%     w(t,:)=exp(logw-max(logw)); 
%     sir_llh=sir_llh+log(mean(w(t,:)))+max(logw(:,1:end));
%     w(t,:)=w(t,:)/sum(w(t,:)); 
% end

% global Dx Du T xx
% Q=tausq*eye(Du); % process noise variance
% Qc=chol(Q); % process noise stdev
% logweights=@(y,x,s)-0.5*(log(2*pi)+log(s)+(y-x).^2./s); % Gaussian log likelihood pdf
% particles=zeros(Dx,N,T); % preallocate
% w=zeros(T,N); % preallocate
% indx=ones(T,N); % preallocate
% particles(1,1,:)=ctraj; % initial value of static particle
% particles(:,2:N,1)=bsxfun(@plus,mu,chol(tausq/(1-phi^2))*randn(Dx,N-1)); % particle initialisation
% logw=logweights(y(1),xx(1,:)*beta,exp(particles(1,:,1))); % log p(y|x)
% w(1,:)=exp(logw-max(logw)); % weights initialisation
% sir_llh=log(mean(w(1,:)))+max(logw); % log likelihood
% w(1,:)=w(1,:)./sum(w(1,:)); % normalisation
% csmchat(:,1)=sum(bsxfun(@times,w(1,:),particles(:,:,1)),2); % weighted mean
% for t=2:T,
%     indx(t,:)=rs_systematic(w(t-1,:),true); % resampling indices
%     particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+Qc*randn(Du,N-1); % state transition
%     logw=logweights(y(t),xx(t,:)*beta,exp(particles(1,:,t))); % log p(y|x)
%     w(t,:)=exp(logw-max(logw)); % update weights
%     sir_llh=sir_llh+log(mean(w(t,:)))+max(logw); % log likelihood
%     w(t,:)=w(t,:)/sum(w(t,:)); % normalisation
%     csmchat(:,t)=sum(bsxfun(@times,w(t,:),particles(:,:,t)),2); % weighted mean
% end