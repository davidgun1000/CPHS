%gen data factor model

phi_true=0.98;
tau_true=0.10225;
rho_true=-0.20;
mu_true=-0.4239;
num_beta=50;
beta_true=0.1*randn(num_beta,1);
T=6000;
%error_term=mvnrnd([0,0],[1,rho_true;rho_true,1],T);
g(1,1)=sqrt(tau_true/(1-phi_true^2))*randn+mu_true;
z(:,1)=mvnrnd(zeros(1,num_beta),eye(num_beta))'; 
y(:,1)=beta_true'*z(:,1)+exp(g(1,1)/2).*randn;
for t=2:T
    g(1,t)=mu_true+phi_true*(g(1,t-1)-mu_true)+rho_true*sqrt(tau_true)*exp(-g(1,t-1)./2).*(y(1,t-1)-beta_true'*z(:,t-1))+sqrt(tau_true)*sqrt(1-rho_true^2)*randn;
    z(:,t)=mvnrnd(zeros(1,num_beta),eye(num_beta))';
    y(:,t)=beta_true'*z(:,t)+exp(g(1,t)./2).*randn;
end

%y=exp(g(:,1)./2).*error_term(:,2);
ctraj_true=g;
save('leverage_6000_cov.mat','y','T','ctraj_true','phi_true','tau_true','rho_true','mu_true','beta_true','z');












%sigma=diag([sigma1,sigma2,sigma3]);
%f=exp(g(:,1)./2).*randn(T,1);
%epsilon=mvnrnd(zeros(1,3),sigma,T);
%epsilon=epsilon';
%y=B*f'+epsilon';
%save('factor_sim_data_mu.mat','y','T');

