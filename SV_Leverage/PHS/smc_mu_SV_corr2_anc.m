%function [particles,w,indx,sir_llh]=smc_mu_SV_corr(y,phi,tau,mu,N,T,u1,u1_res)
function [particles,w,indx,sir_llh]=smc_mu_SV_corr2_anc(y,phi,tau,mu,N,T,u1,u1_res)
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
       %[indx(t,:)]=resampleMultinomialcorr2_anc(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
       %[indx(t,:)]=rs_multinomial_corr3(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
       [indx(t,:)]=rs_multinomial_corr2(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
       %[indx(t,:)]=rs_multinomial_corr(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
       particles(:,1:N,t)=mu+phi*(particles(:,indx(t,:),t-1)-mu)+sqrt(tau)*u1(:,t)';
       logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)).^2);
       w(t,:)=exp(logw-max(logw));
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end

end

%[indx(t,:)]=rs_systematic_corr(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)',cycle(:,t)');
       

% tau=exp(theta_tau);
% t=1;
% particles=zeros(1,N,T);
% w=zeros(T,N);
% indx=zeros(T,N);
% particles(1,1,:)=ctraj;
% particles(:,2:N,1)=sqrt(tau/(1-phi^2))*randn(1,N-1)+mu;
% logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
% w(1,:)=exp(logw-max(logw));
% sir_llh=log(mean(w(1,:)))+max(logw);
% w(1,:)=w(1,:)./sum(w(1,:));
% 
% for t=2:T
%     indx(t,:)=rs_systematic_cond(w(t-1,:));
%     particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+sqrt(tau)*randn(1,N-1);
%     logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*(y(1,t).^2);
%     w(t,:)=exp(logw-max(logw)); 
%     sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
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