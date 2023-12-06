function sirhat=backward_simulation_SV_leverage_univ(particles,w,indx,phi,tau,mu,rho,y)
% backward simulation algorithm
T=size(w,1);
N=size(w,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1,
    eps=(y(1,t)).*exp(-particles(:,:,t)./2);
    log_weight_backward=log(w(t,:))-0.5*log(2*pi)-0.5*log(tau*(1-rho^2))-0.5.*(1/(tau*(1-rho^2))).*((sirhat(:,t+1)-mu-phi.*(particles(:,:,t)-mu)-rho*sqrt(tau).*eps).^2);
    w_backward=exp(log_weight_backward-max(log_weight_backward));
    w_backward=w_backward./sum(w_backward);
    indx(t,1)=find(rand(1) < cumsum(w_backward),1,'first');
    sirhat(:,t)=particles(:,indx(t,1),t);
end
