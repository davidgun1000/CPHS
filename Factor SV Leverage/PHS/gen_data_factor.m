%gen data factor model

dim_y=28;
num_fact=1;
phi_true=0.98;
tau_true=0.05;
mu_true=0.01;
B_true=[1;0.5*ones(27,1)];
T=3000;
t=1;
h=sqrt(tau_true/(1-phi_true^2))*randn(dim_y,1)+mu_true;
ctraj_idio_true(:,t)=h;

hf=sqrt(tau_true/(1-phi_true^2))*randn(num_fact,1);
ctraj_factor_true(:,t)=hf;
for t=2:T
    h=mu_true+phi_true*(h-mu_true)+sqrt(tau_true)*randn(dim_y,1);
    ctraj_idio_true(:,t)=h;    
    hf=phi_true*hf+sqrt(tau_true)*randn(num_fact,1);
    ctraj_factor_true(:,t)=hf;
end

fact_true=exp(ctraj_factor_true'./2).*randn(T,1);
epsilon=exp(ctraj_idio_true./2).*randn(dim_y,T);
y=B_true*fact_true'+epsilon;

save('datasim_dimy28_T3000.mat','y','fact_true','ctraj_idio_true','ctraj_factor_true','T','B_true','phi_true','tau_true','mu_true');


% h=sqrt(tau/(1-phi^2))*randn+mu;
% g(1,1)=h;
% for t=1:T
%     h=mu+phi*(h-mu)+sqrt(tau)*randn;
%     g(t,1)=h;
% end
% sigma=diag([sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,sigma9,sigma10,...
%     sigma11,sigma12,sigma13,sigma14,sigma15,sigma16,sigma17,sigma18,sigma19,sigma20,...
%     sigma21,sigma22,sigma23,sigma24,sigma25,sigma26,sigma27,sigma28,sigma29,sigma30,...
%     sigma31,sigma32,sigma33,sigma34,sigma35,sigma36,sigma37,sigma38,sigma39,sigma40,...
%     sigma41,sigma42,sigma43,sigma44,sigma45,sigma46,sigma47,sigma48,sigma49,sigma50]);
% f=exp(g(:,1)./2).*randn(T,1);
% epsilon=mvnrnd(zeros(1,50),sigma,T);
% %epsilon=epsilon';
% y=B*f'+epsilon';
% save('factor_sim_data_mu50.mat','y','T');

