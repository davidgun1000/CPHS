function sirhat=backward_simulation_SV(particles,w,indx,phi,tau,mu)
%keyboard
%Backward simulation algorithm
T=size(w,1);
N=size(w,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1,
    log_weight_backward=log(w(t,:))-0.5*log(2*pi)-0.5*log(tau)-0.5.*(1/(tau)).*((sirhat(:,t+1)-mu-phi.*(particles(:,:,t)-mu)).^2);
    w_backward=exp(log_weight_backward-max(log_weight_backward));
    w_backward=w_backward./sum(w_backward);
    indx(t,1)=find(rand(1) < cumsum(w_backward),1,'first');
    sirhat(:,t)=particles(:,indx(t,1),t);
end
