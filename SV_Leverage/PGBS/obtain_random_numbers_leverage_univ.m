function [u1_used]=obtain_random_numbers_leverage_univ(phi,tau,mu,rho,ctraj,N,T,y)

for t=1:T
    if t==1
        u1(1,t)=sqrt((1-phi^2)/tau)*(ctraj(1,1)-mu);
    else
        eps=(y(1,t-1)).*exp(-ctraj(1,t-1)./2);
        u1(1,t)=(ctraj(1,t)-mu-phi*(ctraj(1,t-1)-mu)-rho*sqrt(tau)*eps)/(sqrt(tau)*sqrt(1-rho^2));
    end
end

u1_used=[u1;randn(N-1,T)]; 
end
