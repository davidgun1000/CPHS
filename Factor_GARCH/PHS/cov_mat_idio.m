%covmat_idio
dim_y=20;

for s=1:dim_y
    theta=[log(Post.a(2001:end,s)),log(Post.tau(2001:end,s)),Post.mu(2001:end,s)];
    covmat_idio(:,:,s)=cov(theta);
end
