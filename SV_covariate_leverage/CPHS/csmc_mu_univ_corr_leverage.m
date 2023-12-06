function [particles,w,indx,sir_llh,u1,u_res]=csmc_mu_univ_corr_leverage(y,phi,tau,mu,rho,beta,z,N,T,ctraj)

%conditional sequential Monte Carlo algorithm
u1=zeros(N,T);
t=1;
particles=zeros(1,N,T);
w=zeros(T,N);
indx=zeros(T,N);
particles(1,1,:)=ctraj;
u1(1,t)=sqrt((1-phi^2)/tau)*(ctraj(1,1)-mu);
u1(2:N,t)=randn(N-1,1); 
particles(:,2:N,1)=sqrt(tau/(1-phi^2))*u1(2:end,t)'+mu;
logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-beta'*z(:,t)).^2);
w(t,:)=exp(logw-max(logw));
sir_llh=log(mean(w(t,:)))+max(logw);
w(t,:)=w(t,:)./sum(w(t,:)); 
for t=2:T
    eps_ctraj = ((y(1,t-1)-beta'*z(:,t-1)).*exp(-ctraj(1,t-1)./2));
    u1(1,t)=(ctraj(1,t)-mu-phi*(ctraj(1,t-1)-mu)-rho*sqrt(tau)*eps_ctraj)/(sqrt(tau)*sqrt(1-rho^2));
    u1(2:N,t)=randn(N-1,1);
    u=sort(rand(1,N-1));    
    [indx(t,:),u_res(:,t)]=rs_multinomial_cond_corr2(particles(:,1:N,t-1),w(t-1,:),u);
    eps_t=(y(1,t-1)-beta'*z(:,t-1)).*exp(-particles(:,indx(t,2:N),t-1)./2);
    particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+rho*sqrt(tau).*eps_t+sqrt(1-rho^2).*sqrt(tau).*u1(2:end,t)';
    logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-beta'*z(:,t)).^2);
    w(t,:)=exp(logw-max(logw)); 
    sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:));
end
end



% t=1;
%    particles=zeros(1,N,T);
%    w=zeros(T,N);
%    indx=zeros(T,N);
%    particles(1,1,:)=ctraj;
%    particles(:,2:N,1)=sqrt(tau/(1-phi^2))*u1(2:end,t)'+mu;
%    logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-beta'*z(:,t)).^2);
%    w(t,:)=exp(logw-max(logw));
%    sir_llh=log(mean(w(t,:)))+max(logw);
%    w(t,:)=w(t,:)./sum(w(t,:)); 
%    for t=2:T
%        u=sort(rand(1,N-1));            
%        [indx(t,:),u_res(:,t)]=rs_multinomial_cond_corr2(particles(:,1:N,t-1),w(t-1,:),u);
%        eps_t=(y(1,t-1)-beta'*z(:,t-1)).*exp(-particles(:,indx(t,2:N),t-1)./2);
%        particles(:,2:N,t)=mu+phi*(particles(:,indx(t,2:N),t-1)-mu)+rho*sqrt(tau).*eps_t+sqrt(1-rho^2).*sqrt(tau).*u1(2:end,t)';
%        logw=-0.5*log(2*pi)-0.5.*log(exp(particles(1,:,t)))-0.5.*(1./exp(particles(1,:,t))).*((y(1,t)-beta'*z(:,t)).^2);
%        w(t,:)=exp(logw-max(logw)); 
%        sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
%        w(t,:)=w(t,:)./sum(w(t,:));
%    end