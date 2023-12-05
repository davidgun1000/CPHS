function [sirhat]=backward_simulation_GARCH_approx(particles,w,indx,a,tau,mu)
T=size(w,1);
N=size(w,2);
h=1/10;
M=10;
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1,
    if mod(t,M)==1
       %log_weight_backward=log(w(t,:))+(-0.5*log(2*pi)-0.5*log(tau*h)-0.5.*(1/(tau*h)).*((sirhat(:,t+1)-particles(:,:,t)-h.*a.*(mu-particles(:,:,t))).^2));            
       %z=z+h.*(a.*(mu-exp(z)).*exp(-z)-(tau/2))+sqrt(tau)*sqrt(h)*randn(1,N-1);
       
       log_weight_backward=log(w(t,:))+(-0.5*log(2*pi)-0.5*log(tau*h)-0.5.*(1/(tau*h)).*((sirhat(:,t+1)-particles(:,:,t)-h.*(a.*(mu-exp(particles(:,:,t))).*(exp(-particles(:,:,t)))-tau/2)).^2));
       w_backward=exp(log_weight_backward-max(log_weight_backward));
       w_backward=w_backward./sum(w_backward);
       indx_choose=find(rand(1) < cumsum(w_backward),1,'first');
       sirhat(:,t)=particles(:,indx_choose,t);
       outndx(t)=indx(t,indx_choose);
    else
       outndx(t)=indx(t+1,outndx(t+1));
       sirhat(:,t)=particles(:,outndx(t),t);
    end
    
end


end
