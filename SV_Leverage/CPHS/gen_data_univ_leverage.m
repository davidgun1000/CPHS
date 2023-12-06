%gen data factor model

phi=0.98;
tau=0.10225;
rho=-0.4592;
mu=-0.4239;
T=6000;
error_term=mvnrnd([0,0],[1,rho;rho,1],T);
h=sqrt(tau/(1-phi^2))*error_term(1,1)+mu;
g(1,1)=h;
for t=2:T
    h=mu+phi*(h-mu)+sqrt(tau)*error_term(t,1);
    g(t,1)=h;
end

y=exp(g(:,1)./2).*error_term(:,2);
ctraj_true=g;
save('leverage_6000.mat','y','T','ctraj_true');












%sigma=diag([sigma1,sigma2,sigma3]);
%f=exp(g(:,1)./2).*randn(T,1);
%epsilon=mvnrnd(zeros(1,3),sigma,T);
%epsilon=epsilon';
%y=B*f'+epsilon';
%save('factor_sim_data_mu.mat','y','T');

