%IACT_mean

IACT_B1=sum(IACT_used(Post.B1(1001:end,:)));
IACT_tau=sum(IACT_used(Post.tau(1001:end,:)));
IACT_a=sum(IACT_used(Post.a(1001:end,:)));
IACT_mu=sum(IACT_used(Post.mu(1001:end,:)));
IACT_phi=sum(IACT_used(Post.phi(1001:end,:)));
IACT_tau_factor=sum(IACT_used(Post.tau_factor(1001:end,:)));
IACT_mean=(IACT_B1+IACT_tau+IACT_a+IACT_mu+IACT_phi+IACT_tau_factor)/82;