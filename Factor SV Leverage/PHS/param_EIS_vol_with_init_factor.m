function [estimator]=param_EIS_vol_with_init_factor(y,B,ft,phi,tau,mu,N,T,rand_number_CRN,est_init,num,dim_y)

R=size(rand_number_CRN,2);
num_EIS=4;

step=[1;1;1;1;1;1;1;1;1;1];

if num<=dim_y
   var_t=1/(((1-phi^2)/tau)+est_init(1,3));
   mean_t=var_t*(est_init(1,2)+((1-phi^2)/tau)*mu);
   particles(1,:)=sqrt(var_t)*rand_number_CRN(1,:)+mean_t;
   for t=2:T
       var_t=1/((1/tau)+est_init(t,3));
       mean_t=var_t.*(est_init(t,2)+(1/tau).*(mu+phi.*(particles(t-1,:)-mu))');
       particles(t,:)=sqrt(var_t)*rand_number_CRN(t,:)+mean_t';
   end
    
   log_w_1=step(1,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(T,:)))-0.5.*(1./exp(particles(T,:))).*((y(1,T)-B(1,:)*ft(:,T)).^2));
   lw=log_w_1;
   lw=lw';
   X_w=[ones(R,1),particles(T,:)',-0.5.*((particles(T,:)').^2)];
   if sum(isinf(lw))>0 | sum(isnan(lw))>0
     id=isinf(lw) | isnan(lw);
     id=1-id;
     id=logical(id);
     lw=lw(id,1);
     X_w=X_w(id,:); 
  end
  [est]=regress(lw,X_w); 
  [est]=real(est);
  %estimator(T,:)=est';
  estimator(T,:)=est_init(T,:)+0.8.*(est'-est_init(T,:));
  iter=T-1;
  
  while iter>=1
       log_w_1=step(1,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(iter,:)))-0.5.*(1./exp(particles(iter,:))).*((y(1,iter)-B(1,:)*ft(:,iter)).^2));
       var_t=1/((1/tau)+estimator(iter+1,3));
       mean_t=var_t.*(estimator(iter+1,2)+(1/tau).*(mu+phi.*(particles(iter,:)-mu))');
       log_int_const1=0.5*log(tau/var_t);
       log_int_const2=(0.5).*(1/tau).*((mu+phi.*(particles(iter,:)-mu))').^2;
       log_int_const3=-0.5.*(1/var_t).*(mean_t.^2);
       log_int_const=log_int_const1+log_int_const2+log_int_const3;
       lw=log_w_1'-log_int_const;
       X_w=[ones(R,1),particles(iter,:)',-0.5.*((particles(iter,:)').^2)];
       if sum(isinf(lw))>0 | sum(isnan(lw))>0
          id=isinf(lw) | isnan(lw);
          id=1-id;
          id=logical(id);
          lw=lw(id,1);
          X_w=X_w(id,:); 
       end
       [est]=regress(lw,X_w);
       [est]=real(est);
       %estimator(iter,:)=est';
       estimator(iter,:)=est_init(iter,:)+0.8.*(est'-est_init(iter,:));
       iter=iter-1;       
  end
    
  for i=2:num_EIS
      var_t=1/(((1-phi^2)/tau)+estimator(1,3));
      mean_t=var_t*(estimator(1,2)+((1-phi^2)/tau)*mu);
      particles(1,:)=sqrt(var_t)*rand_number_CRN(1,:)+mean_t;
      for t=2:T
         var_t=1/((1/tau)+estimator(t,3));
         mean_t=var_t.*(estimator(t,2)+(1/tau).*(mu+phi.*(particles(t-1,:)-mu))');
         particles(t,:)=sqrt(var_t)*rand_number_CRN(t,:)+mean_t';
      end
      
      log_w_1=step(i,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(T,:)))-0.5.*(1./exp(particles(T,:))).*((y(1,T)-B(1,:)*ft(:,T)).^2));
      lw=log_w_1;
      lw=lw';
      X_w=[ones(R,1),particles(T,:)',-0.5.*((particles(T,:)').^2)];
      if sum(isinf(lw))>0 | sum(isnan(lw))>0
         id=isinf(lw) | isnan(lw);
         id=1-id;
         id=logical(id);
         lw=lw(id,1);
         X_w=X_w(id,:); 
     end
     [est]=regress(lw,X_w);
     %[est]=robustfit(X_w,lw);
     est=real(est);
     estimator(T,:)=estimator(T,:)+0.8.*(est'-estimator(T,:));
     %estimator(T,:)=est';
     iter=T-1; 
      
     while iter>=1
           log_w_1=step(i,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(iter,:)))-0.5.*(1./exp(particles(iter,:))).*((y(1,iter)-B(1,:)*ft(:,iter)).^2));
           var_t=1/((1/tau)+estimator(iter+1,3));
           mean_t=var_t.*(estimator(iter+1,2)+(1/tau).*(mu+phi.*(particles(iter,:)-mu))');
           log_int_const1=0.5*log(tau/var_t);
           log_int_const2=0.5.*(1/tau).*((mu+phi.*(particles(iter,:)-mu))').^2;
           log_int_const3=-0.5.*(1/var_t).*(mean_t.^2);
           log_int_const=log_int_const1+log_int_const2+log_int_const3;
           lw=log_w_1'-log_int_const;
           X_w=[ones(R,1),particles(iter,:)',-0.5.*((particles(iter,:)').^2)];
           if sum(isinf(lw))>0 | sum(isnan(lw))>0
            id=isinf(lw) | isnan(lw);
            id=1-id;
            id=logical(id);
            lw=lw(id,1);
            X_w=X_w(id,:); 
           end 
           [est]=regress(lw,X_w);
           est=real(est);
           %estimator(iter,:)=est';
           estimator(iter,:)=estimator(iter,:)+0.8.*(est'-estimator(iter,:));
           iter=iter-1;
     end
      
      
  end
  estimator=estimator(:,1:3);  

else
    var_t=1/(((1-phi^2)/tau)+est_init(1,3));
    mean_t=var_t*(est_init(1,2)+((1-phi^2)/tau)*mu);
    particles(1,:)=sqrt(var_t)*rand_number_CRN(1,:)+mean_t;
    for t=2:T
        var_t=1/((1/tau)+est_init(t,3));
        mean_t=var_t.*(est_init(t,2)+(1/tau).*(mu+phi.*(particles(t-1,:)-mu))');
        particles(t,:)=sqrt(var_t)*rand_number_CRN(t,:)+mean_t';
    end
    log_w_1=step(1,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(T,:)))-((ft(1,T).^2)./(2.*exp(particles(T,:)))));
    lw=log_w_1;
    lw=lw';
    X_w=[ones(R,1),particles(T,:)',-0.5.*((particles(T,:)').^2)];
    if sum(isinf(lw))>0 | sum(isnan(lw))>0
     id=isinf(lw) | isnan(lw);
     id=1-id;
     id=logical(id);
     lw=lw(id,1);
     X_w=X_w(id,:); 
    end
    [est]=regress(lw,X_w);
    [est]=real(est);
    %estimator(T,:)=est';
    estimator(T,:)=est_init(T,:)+0.8.*(est'-est_init(T,:));
    iter=T-1;
    
    while iter>=1
         log_w_1=step(1,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(iter,:)))-((ft(1,iter).^2)./(2.*exp(particles(iter,:)))));
         var_t=1/((1/tau)+estimator(iter+1,3));
         mean_t=var_t.*(estimator(iter+1,2)+(1/tau).*(mu+phi.*(particles(iter,:)-mu))');
         log_int_const1=0.5*log(tau/var_t);
         log_int_const2=(0.5).*(1/tau).*((mu+phi.*(particles(iter,:)-mu))').^2;
         log_int_const3=-0.5.*(1/var_t).*(mean_t.^2);
         log_int_const=log_int_const1+log_int_const2+log_int_const3;
         lw=log_w_1'-log_int_const;
         X_w=[ones(R,1),particles(iter,:)',-0.5.*((particles(iter,:)').^2)];
         if sum(isinf(lw))>0 | sum(isnan(lw))>0
            id=isinf(lw) | isnan(lw);
            id=1-id;
            id=logical(id);
            lw=lw(id,1);
            X_w=X_w(id,:); 
         end
         [est]=regress(lw,X_w);
         [est]=real(est);
         %estimator(iter,:)=est';
         estimator(iter,:)=est_init(iter,:)+0.8.*(est'-est_init(iter,:));
         iter=iter-1;
    
    end
    
    for i=2:num_EIS
        var_t=1/(((1-phi^2)/tau)+estimator(1,3));
        mean_t=var_t*(estimator(1,2)+((1-phi^2)/tau)*mu);
        particles(1,:)=sqrt(var_t)*rand_number_CRN(1,:)+mean_t;
    for t=2:T
        var_t=1/((1/tau)+estimator(t,3));
        mean_t=var_t.*(estimator(t,2)+(1/tau).*(mu+phi.*(particles(t-1,:)-mu))');
        particles(t,:)=sqrt(var_t)*rand_number_CRN(t,:)+mean_t';
    end
        
    log_w_1=step(i,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(T,:)))-((ft(1,T).^2)./(2.*exp(particles(T,:)))));    
    lw=log_w_1;
    lw=lw';
    X_w=[ones(R,1),particles(T,:)',-0.5.*((particles(T,:)').^2)];     
    if sum(isinf(lw))>0 | sum(isnan(lw))>0
    id=isinf(lw) | isnan(lw);
    id=1-id;
    id=logical(id);
    lw=lw(id,1);
    X_w=X_w(id,:); 
    end    
    [est]=regress(lw,X_w);    
    est=real(est);
    %estimator(T,:)=est';
    estimator(T,:)=estimator(T,:)+0.8.*(est'-estimator(T,:));
    iter=T-1;    
        
    while iter>=1
         log_w_1=step(i,1)*(-0.5*log(2*pi)-0.5*log(exp(particles(iter,:)))-((ft(1,iter).^2)./(2.*exp(particles(iter,:)))));  
         var_t=1/((1/tau)+estimator(iter+1,3));
         mean_t=var_t.*(estimator(iter+1,2)+(1/tau).*(mu+phi.*(particles(iter,:)-mu))');
         log_int_const1=0.5*log(tau/var_t);
         log_int_const2=0.5.*(1/tau).*((mu+phi.*(particles(iter,:)-mu))').^2;
         log_int_const3=-0.5.*(1/var_t).*(mean_t.^2);
         log_int_const=log_int_const1+log_int_const2+log_int_const3;
         lw=log_w_1'-log_int_const;
         X_w=[ones(R,1),particles(iter,:)',-0.5.*((particles(iter,:)').^2)];
         if sum(isinf(lw))>0 | sum(isnan(lw))>0
            id=isinf(lw) | isnan(lw);
            id=1-id;
            id=logical(id);
            lw=lw(id,1);
            X_w=X_w(id,:); 
         end
         [est]=regress(lw,X_w);
         est=real(est);
         %estimator(iter,:)=est';
         estimator(iter,:)=estimator(iter,:)+0.8.*(est'-estimator(iter,:));
         iter=iter-1;
    end
    end
    estimator=estimator(:,1:3);
end

end