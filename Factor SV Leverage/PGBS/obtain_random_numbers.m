function [u1_used]=obtain_random_numbers(phi,tau,mu,ctraj,N,T)

for t=1:T
    if t==1
        u1(1,t)=sqrt((1-phi^2)/tau)*(ctraj(1,1)-mu);
    else
        u1(1,t)=(ctraj(1,t)-mu-phi*(ctraj(1,t-1)-mu))/sqrt(tau);
    end
end

u1_used=[u1;randn(N-1,T)]; 
end
